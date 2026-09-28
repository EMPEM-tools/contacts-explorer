# Contacts Explorer

Explore inter-chain contacts in a protein structure, entirely in your browser.

Upload a `.pdb` file, pick two chains, and get their interface: van der Waals
contacts, hydrogen bonds, salt bridges, π-cations, π/T-stacking and
hydrophobic contacts — with a live 3D view, sortable tables, CSV download,
and ready-to-paste ChimeraX/PyMOL residue selections.

**Nothing is uploaded anywhere.** Your file is read locally with the
browser's FileReader API and all calculations run client-side in JavaScript.
No server, no tracking, no data leaves your machine.

Interaction criteria follow the geometric definitions documented by
[getcontacts](https://getcontacts.github.io/interactions.html)
(sb 4.0 Å · pc 6.0 Å/60° · ps 7.0 Å/30°/45° · ts 5.0 Å/30°/45° · hb 3.5 Å ·
vdw/hp summed radii + 0.5 Å), reimplemented independently in JavaScript.

## Use

Open [the GitHub Page](.) and drop in a PDB. That's it.

## Local development

Static site — open `index.html` directly or serve the folder:
`python3 -m http.server`
