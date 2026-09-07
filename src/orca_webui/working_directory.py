"""Working-directory column: directory selection, file browser, and file viewers.

Owns the left-hand column of the UI. A "working directory" is a subfolder under
``data/`` in which all structure files, ORCA inputs/outputs, and exports live. This
module also holds the two shared ``gr.State`` objects (current path and file list)
that the other tabs subscribe to; any handler that changes files on disk returns a
fresh file list into that state so every tab's dropdowns stay in sync.
"""
import os
import shutil
import tempfile
import time
import uuid
from functools import partial
from pathlib import Path
import pandas as pd
import gradio as gr
import nglview
from .utils import get_files_in_working_directory, mol_from_structure_file, resolve_working_file

def get_working_directories(base_path="./data/") -> list[str]:
    """List the existing working-directory names (immediate subfolders of ``data/``)."""
    base = Path(base_path).resolve()
    if not base.is_dir():
        return []
    return sorted((d.name for d in base.iterdir()
                   if d.is_dir() and d.resolve().parent == base), key=str.casefold)

def on_open_working_directory(working_directory, base_path="./data/"):
    """Create (if needed) and open ``working_directory`` under ``data/``.

    Returns updates for the directory dropdown, the path/file-list states, and the
    two buttons that become interactive once a directory is open. Warns and no-ops
    when the name is blank.
    """
    try:
        if not isinstance(working_directory, str) or not working_directory.strip():
            raise ValueError("Please specify a working directory.")
        working_directory = working_directory.strip()
        if working_directory in {".", ".."} or any(c in working_directory for c in "/\\\0"):
            raise ValueError("Use a directory name without path separators.")
        base = Path(base_path).resolve()
        directory = (base / working_directory).resolve()
        if directory.parent != base:
            raise ValueError("The working directory must be inside the data directory.")
        directory.mkdir(parents=True, exist_ok=True)
        files = get_files_in_working_directory(str(directory))
        return gr.update(choices=get_working_directories(base), value=working_directory), str(directory), files, gr.update(interactive=True), gr.update(interactive=True)
    except (OSError, ValueError) as exc:
        gr.Warning(str(exc))
        # Preserve the active directory and its controls when opening another fails.
        return tuple(gr.update() for _ in range(5))

def on_file_list_change(working_directory_path) -> pd.DataFrame:
    """Build the file-browser DataFrame (File / Type / Modified) for the directory.

    Each file is classified by extension into a human-readable type, and rows are
    ordered by actual modification time (newest first). Yields an empty table when no
    working directory is open, since this fires on the shared state's ``.change`` and
    that state is ``None`` before the first directory is opened.
    """
    if not working_directory_path:
        return pd.DataFrame(columns=["File", "Type", "Modified"])

    files = get_files_in_working_directory(working_directory_path)

    # Update the file dataframe
    file_info = []
    for f in files:
        if f.endswith('.xyz') or f.endswith('.pdb') or f.endswith('.mol') or f.endswith('.mol2'):
            file_type = "Structure file"
        elif f.endswith('.inp'):
            file_type = "Input file"
        elif f.endswith('.gbw'):
            file_type = "General binary wavefunction file"
        elif f.endswith('.scfp'):
            file_type = "SCF density matrix file"
        elif f.endswith('.densities') or f.endswith('.densitiesinfo'):
            file_type = "Density data file"
        elif f.endswith('.engrad'):
            file_type = "Energy and gradients file"
        elif f.endswith('.log'):
            file_type = "Log file"
        elif f.endswith('.cube'):
            file_type = "Cube file"
        elif f.endswith('.txt'):
            file_type = "Text file"
        else:
            file_type = "Other File"
        
        file_path = os.path.join(working_directory_path, f)
        try:
            mtime = os.path.getmtime(file_path)
        except OSError:
            continue  # The file may disappear while ORCA is cleaning up scratch files.
        file_info.append([f, file_type, mtime])
    # Sort by the numeric modification time (newest first), then format for display.
    # Sorting the ctime *string* would order by weekday name, not chronologically.
    file_info.sort(key=lambda x: x[2], reverse=True)
    file_info = [[f, file_type, time.ctime(mtime)] for f, file_type, mtime in file_info]
    file_df = pd.DataFrame(file_info, columns=["File", "Type", "Modified"])

    return file_df

def on_select_file(evt: gr.SelectData):
    """Handle a row click in the file browser.

    Returns the selected file name plus two derived selections (a viewable structure
    file and a viewable text file, each ``None`` when the type doesn't apply) and
    enables the delete button. These drive which viewer buttons become active.
    """
    selected_file_name = evt.row_value[0]
    if selected_file_name.endswith('.xyz') or selected_file_name.endswith('.pdb') or selected_file_name.endswith('.mol') or selected_file_name.endswith('.mol2') or selected_file_name.endswith('.log'):
        selected_structure_file = selected_file_name
    else:
        selected_structure_file = None
    if selected_file_name.endswith('.xyz') or selected_file_name.endswith('.pdb') or selected_file_name.endswith('.mol') or selected_file_name.endswith('.mol2') or selected_file_name.endswith('.inp') or selected_file_name.endswith('.log') or selected_file_name.endswith('.engrad') or selected_file_name.endswith('.txt'):
        selected_text_file = selected_file_name
    else:
        selected_text_file = None
    
    return selected_file_name, selected_structure_file, selected_text_file, gr.update(interactive=True)

def on_selected_structure_file_state_change(state):
    """Enable the "View Structure File" button only when a structure file is selected."""
    return gr.update(interactive=(state is not None))

def on_selected_text_file_state_change(state):
    """Enable the "View Text File" button only when a text-viewable file is selected."""
    return gr.update(interactive=(state is not None))

def on_upload_file(working_directory_path, file_path):
    """Copy an uploaded file into the working directory; return the refreshed file list."""
    try:
        if not file_path:
            raise ValueError("Please select a file to upload.")
        destination = resolve_working_file(working_directory_path, os.path.basename(file_path))
        if os.path.lexists(destination):
            raise ValueError("A file with this name already exists. Rename the upload first.")
        shutil.copy2(file_path, destination)
    except (OSError, TypeError, ValueError) as exc:
        gr.Warning(f"Error uploading file!\n{exc}")
    return get_files_in_working_directory(working_directory_path)

def on_delete_file(working_directory_path, selected_file_name):
    """Delete the selected file; return the refreshed file list (warns on error)."""
    if selected_file_name is None:
        return get_files_in_working_directory(working_directory_path)
    
    try:
        file_path = resolve_working_file(working_directory_path, selected_file_name)
        os.remove(file_path)
        status = "File deleted successfully."
    except Exception as exc:
        status = "Error deleting file!\n" + str(exc)
    gr.Warning(status)
    
    return get_files_in_working_directory(working_directory_path)

def on_clean_working_directory(working_directory_path):
    """Remove scratch/leftover files (ORCA temp, ``.bibtex``, ``Zone.Identifier``, etc.).

    Returns the refreshed file list; warns on error.
    """
    try:
        if not working_directory_path or not os.path.isdir(working_directory_path):
            raise ValueError("Please open a working directory first.")
        files_to_clean = [f for f in os.listdir(working_directory_path) if f.startswith('#') or f.endswith(".tmp") or f.endswith(".bibtex") or f.endswith("Zone.Identifier") or ".tmp." in f]
        for f in files_to_clean:
            file_path = resolve_working_file(working_directory_path, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
        status = "Working directory cleaned successfully."
    except Exception as exc:
        status = "Error cleaning working directory!\n" + str(exc)
    gr.Warning(status)
    
    return get_files_in_working_directory(working_directory_path)

def on_view_structure_file(working_directory_path, file_name, static_dir="./static"):
    """Render a 3D viewer for the selected structure/output file.

    Loads the geometry with :func:`mol_from_structure_file` — the same loader the
    Calculation tab uses to build ORCA inputs, so the atom indices labelled here are the
    indices of the generated ``.inp`` coordinate block. Builds an nglview widget, writes
    it to a unique file under ``static/``, and returns an ``<iframe>`` pointing at
    it. Returns ``None`` and warns on error.
    """
    try:
        file_path = resolve_working_file(working_directory_path, file_name)
        mol = mol_from_structure_file(file_path)

        # Create the NGL view widget
        view = nglview.show_rdkit(mol)

        # Label each atom with its position in the file / generated .inp coordinate
        # block, 0-based, e.g. ethanol -> C0 C1 O2 H3..H8.
        #
        # labelType="text" with an explicit labelText is what makes this work.
        # labelType="atomname" reads the names out of the PDB block that show_rdkit
        # generates, and RDKit numbers those per element (C1, C2, O1, ...) -- so the
        # third atom of ethanol came out as "O1", the first oxygen, rather than "O2".
        # labelText takes one entry per atom of the selection, in atom order.
        atom_labels = [f"{atom.GetSymbol()}{atom.GetIdx()}" for atom in mol.GetAtoms()]
        view.add_representation("label", labelType="text", labelText=atom_labels,
                                color="black", showBackground=False)

        # Write the widget to HTML
        static_path = Path(static_dir)
        static_path.mkdir(parents=True, exist_ok=True)
        viewer_name = f"structure_{uuid.uuid4().hex}.html"
        viewer_path = resolve_working_file(str(static_path), viewer_name)
        nglview.write_html(viewer_path, [view])

        html = f'<iframe src="/static/{viewer_name}" height="600" width="600" title="NGL View"></iframe>'

        return html
    except Exception as exc:
        gr.Warning("Error!\n" + str(exc))
        return None

def on_view_text_file(working_directory_path, text_file_name):
    """Load a text-viewable file into the editable text viewer.

    Returns updates for the text area (label + content, made editable), save button,
    and the absolute path of the loaded file. Warns and no-ops on read error.
    """
    try:
        text_file_path = resolve_working_file(working_directory_path, text_file_name)
        with open(text_file_path, 'r') as file:
            content = file.read()
        return gr.update(label=f"Text File Viewer - {text_file_name}", value=content, interactive=True), gr.update(interactive=True), text_file_path
    except Exception as exc:
        gr.Warning("Error!\n" + str(exc))
        return gr.update(), gr.update(), gr.update()

def on_save_viewed_text_file(working_directory_path, viewed_file_path, text_content):
    """Save to the file loaded into the editor, even when a different row is selected."""
    try:
        if not viewed_file_path:
            raise ValueError("Please view a text file before saving.")
        name = os.path.basename(viewed_file_path)
        current_path = resolve_working_file(working_directory_path, name)
        if current_path != viewed_file_path or not os.path.isfile(current_path):
            raise ValueError("The viewed file is no longer in the active working directory. Open it again.")
        return on_save_text_file(working_directory_path, name, text_content)
    except (OSError, TypeError, ValueError) as exc:
        gr.Warning(str(exc))
        return get_files_in_working_directory(working_directory_path)

def on_save_text_file(working_directory_path, text_file_name, text_content):
    """Write the text-viewer contents back to the selected file.

    Returns the refreshed file list; warns if no file is selected or on write error.
    """
    if text_file_name is None:
        gr.Warning("Please select a text file to save.")
        return get_files_in_working_directory(working_directory_path)
    
    temporary_path = None
    try:
        text_file_path = resolve_working_file(working_directory_path, text_file_name)
        with tempfile.NamedTemporaryFile(mode='w', dir=working_directory_path, delete=False) as file:
            temporary_path = file.name
            file.write(text_content)
        os.replace(temporary_path, text_file_path)
        status = "File saved successfully."
    except Exception as exc:
        status = "Error saving file!\n" + str(exc)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            try:
                os.remove(temporary_path)
            except OSError:
                pass
    gr.Warning(status)
    
    return get_files_in_working_directory(working_directory_path)

def on_working_directory_change():
    """Discard file selections and editor contents when switching directories."""
    return (None, None, None, None, gr.update(interactive=False), None,
            gr.update(label="Text File Viewer", value="", interactive=False),
            gr.update(interactive=False))

def working_directory_blocks(base_path="./data", static_dir="./static"):
    """Build the working-directory column and wire its events.

    Returns the two shared ``gr.State`` objects — ``(working_directory_path_state,
    working_directory_file_list_state)`` — which the calculation and result tabs
    subscribe to.
    """
    with gr.Column(scale=1):
        working_directory_dropdown = gr.Dropdown(label="Working Directory", choices=get_working_directories(base_path), value="wd", allow_custom_value=True)
        working_directory_path_state = gr.State()
        open_working_directory_button = gr.Button(value="Create/Open Working Directory")
        working_directory_file_list_state = gr.State()
        working_directory_file_dataframe = gr.Dataframe(label="Files in Working Directory", headers=["File", "Type", "Modified"], max_height=360, wrap=True, interactive=False)
        selected_file_state = gr.State()
        selected_structure_file_state = gr.State()
        selected_text_file_state = gr.State()
        viewed_text_file_state = gr.State()
        with gr.Row():
            add_file_upload_button = gr.UploadButton(label="Add File", file_types=[".xyz", ".pdb", ".mol", ".mol2"], interactive=False)
            delete_file_button = gr.Button(value="Delete Selected File", interactive=False)
            clean_working_directory_button = gr.Button(value="Clean Working Directory", interactive=False)
        view_structure_button = gr.Button(value="View Structure File", interactive=False)
        structure_viewer_html = gr.HTML()
        view_text_file_button = gr.Button(value="View Text File", interactive=False)
        text_file_viewer_textarea = gr.TextArea(label="Text File Viewer", lines=20, elem_id="textfile_viewer", interactive=False)
        save_text_file_button = gr.Button(value="Save Text File", interactive=False)
    
    open_directory = partial(on_open_working_directory, base_path=base_path)
    working_directory_dropdown.change(open_directory, working_directory_dropdown, [working_directory_dropdown, working_directory_path_state, working_directory_file_list_state, add_file_upload_button, clean_working_directory_button])
    open_working_directory_button.click(open_directory, working_directory_dropdown, [working_directory_dropdown, working_directory_path_state, working_directory_file_list_state, add_file_upload_button, clean_working_directory_button])
    working_directory_path_state.change(on_working_directory_change, [], [selected_file_state, selected_structure_file_state, selected_text_file_state, viewed_text_file_state, delete_file_button, structure_viewer_html, text_file_viewer_textarea, save_text_file_button])
    working_directory_file_list_state.change(on_file_list_change, working_directory_path_state, working_directory_file_dataframe)
    working_directory_file_dataframe.select(on_select_file, [], [selected_file_state, selected_structure_file_state, selected_text_file_state, delete_file_button])
    selected_structure_file_state.change(on_selected_structure_file_state_change, selected_structure_file_state, view_structure_button)
    selected_text_file_state.change(on_selected_text_file_state_change, selected_text_file_state, view_text_file_button)
    add_file_upload_button.upload(on_upload_file, [working_directory_path_state, add_file_upload_button], working_directory_file_list_state)
    delete_file_button.click(on_delete_file, [working_directory_path_state, selected_file_state], working_directory_file_list_state)
    clean_working_directory_button.click(on_clean_working_directory, working_directory_path_state, working_directory_file_list_state)
    view_structure_button.click(partial(on_view_structure_file, static_dir=static_dir), [working_directory_path_state, selected_structure_file_state], structure_viewer_html)
    view_text_file_button.click(on_view_text_file, [working_directory_path_state, selected_text_file_state], [text_file_viewer_textarea, save_text_file_button, viewed_text_file_state])
    save_text_file_button.click(on_save_viewed_text_file, [working_directory_path_state, viewed_text_file_state, text_file_viewer_textarea], working_directory_file_list_state)

    return working_directory_path_state, working_directory_file_list_state
