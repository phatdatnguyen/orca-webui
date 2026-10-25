"""Offline checks for the packaged editor and safe, atomic asset installation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import stat
import zipfile

import pytest


@pytest.fixture
def vendor(monkeypatch):
    helper = Path(__file__).resolve().parents[1] / "scripts" / "vendor_molecular_editor.py"
    spec = importlib.util.spec_from_file_location("vendor_molecular_editor", helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def no_network(*args, **kwargs):
        pytest.fail("Asset verification and cached installs must not use the network")

    monkeypatch.setattr(module, "urlopen", no_network)
    return module


def test_bundled_editor_is_complete_and_matches_recorded_checksums(vendor):
    assert vendor.installation_is_current(vendor.KETCHER_DIR)
    marker = json.loads((vendor.KETCHER_DIR / vendor.MARKER).read_text())
    actual = {
        path.relative_to(vendor.KETCHER_DIR).as_posix()
        for path in vendor.KETCHER_DIR.rglob("*") if path.is_file()
    }
    assert actual == set(marker["files"]) | {vendor.MARKER}
    for name, expected in vendor.NOTICE_SHA256.items():
        assert vendor.file_sha256(vendor.KETCHER_DIR / name) == expected


@pytest.mark.parametrize("entry_name", [
    "../escape", "nested/../../escape", "/absolute", "C:/escape",
    "C:escape", "nested\\..\\escape", "index.html:alternate-stream",
])
def test_archive_rejects_paths_unsafe_on_posix_or_windows(tmp_path, vendor, entry_name):
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr(entry_name, "must not be extracted")
    staging = tmp_path / "staging"
    staging.mkdir()

    with pytest.raises(ValueError, match="Unsafe editor archive path"):
        vendor.extract_archive(archive, staging)

    assert list(staging.iterdir()) == []
    assert not (tmp_path / "escape").exists()


def test_archive_rejects_symlinks(tmp_path, vendor):
    archive = tmp_path / "symlink.zip"
    link = zipfile.ZipInfo("index.html")
    link.create_system = 3
    link.external_attr = (stat.S_IFLNK | 0o777) << 16
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr(link, "../../outside.html")
    staging = tmp_path / "staging"
    staging.mkdir()

    with pytest.raises(ValueError, match="Unsupported editor archive entry"):
        vendor.extract_archive(archive, staging)

    assert list(staging.iterdir()) == []


def test_checksum_failure_leaves_no_installation_or_staging_files(tmp_path, vendor):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "release.zip").write_bytes(b"not the pinned release")
    destination = tmp_path / "package" / "ketcher"

    with pytest.raises(ValueError, match="Checksum mismatch"):
        vendor.install_ketcher(destination, cache_dir=cache)

    assert not destination.exists()
    assert list(destination.parent.iterdir()) == []


def test_existing_unverified_installation_is_preserved(tmp_path, vendor):
    destination = tmp_path / "ketcher"
    destination.mkdir()
    existing = destination / "index.html"
    existing.write_text("local changes must survive")

    with pytest.raises(RuntimeError, match="no files were changed"):
        vendor.install_ketcher(destination)

    assert list(destination.iterdir()) == [existing]
    assert existing.read_text() == "local changes must survive"


@pytest.fixture
def cached_installation(tmp_path, monkeypatch, vendor):
    """A small release exercises the complete installer without copying 120 MB."""
    cache = tmp_path / "cache"
    cache.mkdir()
    with zipfile.ZipFile(cache / "release.zip", "w") as archive:
        for name in vendor.REQUIRED_FILES:
            if name not in vendor.NOTICE_SHA256:
                archive.writestr(name, f"synthetic asset: {name}")
    monkeypatch.setattr(vendor, "ARCHIVE_SHA256", vendor.file_sha256(cache / "release.zip"))
    checksums = {}
    for name in vendor.NOTICE_SHA256:
        content = f"synthetic notice: {name}".encode()
        (cache / name).write_bytes(content)
        checksums[name] = hashlib.sha256(content).hexdigest()
    monkeypatch.setattr(vendor, "NOTICE_SHA256", checksums)
    destination = tmp_path / "package" / "ketcher"
    assert vendor.install_ketcher(destination, cache_dir=cache) == destination
    return destination


def test_offline_install_includes_notices_and_reuses_verified_assets(vendor, cached_installation):
    assert vendor.installation_is_current(cached_installation)
    assert (cached_installation / "LICENSE").is_file()
    assert (cached_installation / "NOTICE").is_file()
    # A second run needs neither a cache nor downloads and leaves assets in place.
    marker_before = (cached_installation / vendor.MARKER).read_bytes()
    assert vendor.install_ketcher(cached_installation) == cached_installation
    assert (cached_installation / vendor.MARKER).read_bytes() == marker_before


@pytest.mark.parametrize("damage", ["modified", "missing", "outside_symlink"])
def test_verification_rejects_damaged_or_escaping_assets(vendor, cached_installation, tmp_path, damage):
    asset = cached_installation / "index.html"
    original = asset.read_bytes()
    if damage == "modified":
        asset.write_bytes(original + b"changed")
    else:
        asset.unlink()
        if damage == "outside_symlink":
            outside = tmp_path / "outside.html"
            outside.write_bytes(original)
            asset.symlink_to(outside)

    assert not vendor.installation_is_current(cached_installation)
    with pytest.raises(RuntimeError, match="no files were changed"):
        vendor.install_ketcher(cached_installation)
