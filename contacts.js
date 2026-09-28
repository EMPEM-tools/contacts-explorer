/* Client-side contact engine.
 *
 * Implements the geometric criteria documented at
 * https://getcontacts.github.io/interactions.html with the CLI defaults from
 * get_static_contacts.py's argparser:
 *   sb  4.0 A between anion (ASP OD1+OD2, GLU OE1+OE2) and cation (LYS NZ,
 *       ARG NH1+NH2) atom groups
 *   pc  6.0 A cation to aromatic-centroid, ring-normal angle <= 60 deg
 *   ps  7.0 A centroid-centroid, normals <= 30 deg, psi <= 45 deg
 *   ts  5.0 A centroid-centroid, |normals angle - 90| <= 30 deg, psi <= 45 deg
 *   hb  3.5 A donor-acceptor, static structure: D-H-A angle >= 180 deg is
 *       the default, which cannot be evaluated without hydrogens; with H
 *       present we require angle(DHA) <= 40 deg from linear (i.e. >=140 deg)
 *       as VMD's "ideal" hbond criterion; without H, donor-acceptor distance
 *       alone
 *   hp  heavy atoms of hydrophobic residues (C/S element) within
 *       vdw(A)+vdw(B)+0.5
 *   vdw any non-H pair within vdw(A)+vdw(B)+0.5
 *
 * Independent implementation of the published criteria (getcontacts itself
 * is GPL; the geometric rules are chemistry). Output shape mirrors the
 * internal tool's so the UI code ports unchanged.
 */

// VMD/ChimeraX-flavoured vdW radii (A) by element
const VDW = { H: 1.10, C: 1.70, N: 1.55, O: 1.52, S: 1.80, P: 1.80,
              F: 1.47, CL: 1.75, BR: 1.85, I: 1.98, NA: 2.27, MG: 1.73,
              K: 2.75, CA: 2.31, ZN: 1.39, FE: 1.95, MN: 2.05, CU: 1.4 };

const THREE2ONE = { ALA:'A', ARG:'R', ASN:'N', ASP:'D', CYS:'C', GLN:'Q', GLU:'E',
  GLY:'G', HIS:'H', ILE:'I', LEU:'L', LYS:'K', MET:'M', PHE:'F', PRO:'P',
  SER:'S', THR:'T', TRP:'W', TYR:'Y', VAL:'V' };

const HYDROPHOBIC = new Set(['ALA','CYS','PHE','GLY','ILE','LEU','MET','PRO','VAL','TRP']);

// aromatic ring atom groups per residue (getcontacts definitions)
const RINGS = {
  PHE: [ ['CG','CD1','CE1','CE2','CD2'] ],        // centroid over CG+CD1+CE1+CE2+CD2
  TYR: [ ['CG','CD1','CE1','CE2','CD2'] ],
  TRP: [ ['CD2','CE3','CZ3','CH2','CZ2'] ],        // six-membered ring only
};

// ---- PDB parsing ----------------------------------------------------------

function parsePDB(text) {
  const atoms = [];          // {chain, resnum, icode, resname, name, x, y, z, element}
  const seen = new Set();    // (chain, resnum, icode) residue dedup -> chain list
  const chains = {};
  const lines = text.split('\n');
  for (const line of lines) {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
    // fixed-width columns; tolerate short lines
    if (line.length < 54) continue;
    const chain = line[21];
    const resnum = parseInt(line.substring(22, 26), 10);
    if (isNaN(resnum)) continue;
    const icode = (line[26] || ' ').trim();
    const resname = line.substring(17, 20).trim();
    const name = line.substring(12, 16).trim();
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (isNaN(x) || isNaN(y) || isNaN(z)) continue;
    let element = line.length >= 78 ? line.substring(76, 78).trim() : '';
    if (!element) element = name.replace(/[^A-Za-z]/g, '').substring(0, 1).toUpperCase();
    if (element.length > 1) element = element[0] + element.substring(1).toLowerCase();
    const at = { chain, resnum, icode, resname, name, x, y, z, element };
    atoms.push(at);
    const key = chain + '|' + resnum + '|' + icode;
    if (!seen.has(key)) {
      seen.add(key);
        (chains[chain] = chains[chain] || []).push({ resnum, icode, resname });
    }
  }
  // water chains excluded from the selectable list (matches the internal tool)
  for (const c of Object.keys(chains)) {
    const waterFrac = chains[c].filter(r => r.resname === 'HOH' || r.resname === 'WAT').length
                    / Math.max(1, chains[c].length);
    if (waterFrac > 0.5) { chains[c]._water = true; }
  }
  return { atoms, chains };
}

// ---- geometry helpers -----------------------------------------------------

const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
function centroid(atoms) {
  let x = 0, y = 0, z = 0;
  for (const a of atoms) { x += a.x; y += a.y; z += a.z; }
  const n = atoms.length;
  return { x: x / n, y: y / n, z: z / n };
}
function cross(u, v) {
  return { x: u.y*v.z - u.z*v.y, y: u.z*v.x - u.x*v.z, z: u.x*v.y - u.y*v.x };
}
function sub(u, v) { return { x: u.x - v.x, y: u.y - v.y, z: u.z - v.z }; }
function dot(u, v) { return u.x*v.x + u.y*v.y + u.z*v.z; }
function norm(u) { return Math.hypot(u.x, u.y, u.z); }
// angle in degrees between two vectors
function angleDeg(u, v) {
  const c = dot(u, v) / (norm(u) * norm(v) || 1);
  return Math.acos(Math.min(1, Math.max(-1, c))) * 180 / Math.PI;
}

// ---- interaction detection ------------------------------------------------

const ANION = { ASP: ['OD1', 'OD2'], GLU: ['OE1', 'OE2'] };
const CATION = { LYS: ['NZ'], ARG: ['NH1', 'NH2'] };

function ringGroups(res) {
  // all aromatic ring groups (centroids + normals) in a residue
  const out = [];
  const def = RINGS[res.resname];
  if (!def) return out;
  for (const names of def) {
    const ats = res.atoms.filter(a => names.includes(a.name));
    if (ats.length === names.length) {
      const c = centroid(ats);
      const n1 = cross(sub(ats[1], ats[0]), sub(ats[2], ats[0]));
      out.push({ center: c, normal: n1 });
    }
  }
  return out;
}

/* Compute contacts between chains. Returns rows in the internal tool's
 * shape: {atom1, chain1, resname1, resnum1, resnum1_label, itype, dist,
 * atom2, chain2, ...}. */
function computeContacts(pdb, chain1, chain2) {
  const _t0 = Date.now();
  const { atoms } = typeof pdb === 'string' ? parsePDB(pdb) : pdb;
  const a1 = atoms.filter(a => a.chain === chain1);
  const a2 = atoms.filter(a => a.chain === chain2);
  const rows = [];
  const push = (at1, at2, itype, d) => rows.push({
    atom1: at1.name, chain1,
    resname1: at1.resname, resnum1: at1.resnum,
    resnum1_label: at1.resnum + at1.icode,
    itype, dist: Math.round(d * 100) / 100,
    atom2: at2.name, chain2,
    resname2: at2.resname, resnum2: at2.resnum,
    resnum2_label: at2.resnum + at2.icode,
  });

  // index residues by (chain,resnum,icode)
  const resKey = a => a.chain + '|' + a.resnum + '|' + a.icode;
  const residues = new Map();
  for (const a of [...a1, ...a2]) {
    const k = resKey(a);
    if (!residues.has(k)) residues.set(k, { chain: a.chain, resnum: a.resnum, icode: a.icode,
                                            resname: a.resname, atoms: [] });
    residues.get(k).atoms.push(a);
  }

  const heavy1 = a1.filter(a => a.element !== 'H');
  const heavy2 = a2.filter(a => a.element !== 'H');
  if (typeof console !== 'undefined' && console.debug) console.debug('parse/filter', Date.now()-_t0);
  const _t1 = Date.now();

  // --- vdw + hp: heavy-atom pairs in range via uniform grid (max reach
  // is sum of two largest radii + eps ~ 4.1 A; cell 8 A with 1-cell halo) ---
  const CELL = 8.0;
  const grid = new Map();
  const cellKey = (x, y, z) => Math.floor(x / CELL) + ',' + Math.floor(y / CELL) + ',' + Math.floor(z / CELL);
  for (const at of heavy1) {
    const k = cellKey(at.x, at.y, at.z);
    if (!grid.has(k)) grid.set(k, []);
    grid.get(k).push(at);
  }
  for (const at2 of heavy2) {
    const cx = Math.floor(at2.x / CELL), cy = Math.floor(at2.y / CELL), cz = Math.floor(at2.z / CELL);
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const cell = grid.get((cx + dx) + ',' + (cy + dy) + ',' + (cz + dz));
      if (!cell) continue;
      for (const at1 of cell) {
        const d = dist(at1, at2);
        const r = (VDW[at1.element] || 1.7) + (VDW[at2.element] || 1.7) + 0.5;
        if (d > r) continue;
        if (HYDROPHOBIC.has(at1.resname) && HYDROPHOBIC.has(at2.resname)
            && (at1.element === 'C' || at1.element === 'S')
            && (at2.element === 'C' || at2.element === 'S')) {
          push(at1, at2, 'hp', d);
        } else {
          push(at1, at2, 'vdw', d);
        }
      }
    }
  }

  if (typeof console !== 'undefined' && console.debug) console.debug('vdw phase', Date.now()-_t1);
  const _t2 = Date.now();
  // --- hb: donor-acceptor pairs (N/O donors+acceptors, with the angle
  // criterion when hydrogens are present) ---
  const isNO = a => a.element === 'N' || a.element === 'O';
  const hOf = {};   // map parent atom key -> attached H atoms
  // grid heavy atoms once for the covalent attachment lookup
  const covGrid = new Map();
  for (const p of [...heavy1, ...heavy2]) {
    const k = Math.floor(p.x / 4) + ',' + Math.floor(p.y / 4) + ',' + Math.floor(p.z / 4);
    if (!covGrid.has(k)) covGrid.set(k, []);
    covGrid.get(k).push(p);
  }
  const atomKey = a => a.chain + '|' + a.resnum + '|' + a.icode + '|' + a.name;
  for (const a of atoms) if (a.element === 'H') {
    // attach to nearest heavy atom (covalent range ~1.3 A) via the grid
    let best = null, bd = 1.4;
    const cx = Math.floor(a.x / 4), cy = Math.floor(a.y / 4), cz = Math.floor(a.z / 4);
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const cell = covGrid.get((cx + dx) + ',' + (cy + dy) + ',' + (cz + dz));
      if (!cell) continue;
      for (const p of cell) {
        const d = dist(a, p);
        if (d < bd) { bd = d; best = p; }
      }
    }
    if (best) (hOf[atomKey(best)] = hOf[atomKey(best)] || []).push(a);
  }
  const hasH = Object.keys(hOf).length > 0;
  const noGrid = new Map();
  for (const at of heavy1.filter(isNO)) {
    const k = cellKey(at.x, at.y, at.z);
    if (!noGrid.has(k)) noGrid.set(k, []);
    noGrid.get(k).push(at);
  }
  for (const at2 of heavy2.filter(isNO)) {
    const cx = Math.floor(at2.x / CELL), cy = Math.floor(at2.y / CELL), cz = Math.floor(at2.z / CELL);
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const cellA = noGrid.get((cx + dx) + ',' + (cy + dy) + ',' + (cz + dz));
      if (!cellA) continue;
      for (const at1 of cellA) {
      const d = dist(at1, at2);
      if (d > 3.5) continue;
      if (hasH) {
        // need an H attached to at1 or at2 forming a near-linear D-H-A
        let ok = false;
        const key = c => c.chain + '|' + c.resnum + '|' + c.icode + '|' + c.name;
        for (const h of (hOf[key(at1)] || [])) {
          const ang = angleDeg(sub(h, at1), sub(at2, at1));   // H->D vs D->A
          if (ang <= 40) { ok = true; break; }
        }
        for (const h of (hOf[key(at2)] || [])) {
          const ang = angleDeg(sub(h, at2), sub(at1, at2));
          if (ang <= 40) { ok = true; break; }
        }
        if (!ok) continue;
      }
      {
      // side-chain / backbone profile like getcontacts (hbbb/hbsb/hbss)
      const bb = a => a.name === 'N' || a.name === 'O' || a.name.startsWith('C');
      const sc1 = !bb(at1), sc2 = !bb(at2);
      const tag = sc1 && sc2 ? 'hbss' : sc1 ? 'hbsb' : 'hbbb';
      push(at1, at2, tag, d);
      }
    }
  }
  }

  if (typeof console !== 'undefined' && console.debug) console.debug('hb phase', Date.now()-_t2);
  const _t3 = Date.now();
  // --- sb: anion atom-group to cation atom-group ---
  const anionAtoms = heavy1.concat(heavy2).filter(a => ANION[a.resname]);
  const cationAtoms = heavy1.concat(heavy2).filter(a => CATION[a.resname]);
  for (const an of anionAtoms) {
    for (const ca of cationAtoms) {
      if (an.chain === ca.chain) continue;
      if ((an.chain === chain1 ? heavy1 : heavy2).indexOf(an) < 0) continue;
      const d = dist(an, ca);
      if (d <= 4.0) {
        // orientation: an1->an2->ca (ASP/GLU symmetric pair) reported at atom level
        push(an.chain === chain1 ? an : ca, an.chain === chain1 ? ca : an, 'sb', d);
      }
    }
  }

  // --- pc / ps / ts over aromatic rings ---
  const rings1 = [], rings2 = [];
  for (const [k, res] of residues) {
    for (const g of ringGroups(res)) (res.chain === chain1 ? rings1 : rings2).push(g);
  }
  const cationList = [...a1, ...a2].filter(a => CATION[a.resname]);
  // pc: cation near a ring in the OTHER chain
  for (const g of [...rings1, ...rings2]) {
    const ringChain = rings1.includes(g) ? chain1 : chain2;
    for (const ca of cationList) {
      if (ca.chain === ringChain) continue;
      const v = sub({ x: ca.x, y: ca.y, z: ca.z }, g.center);
      const d = norm(v);
      if (d > 6.0) continue;
      if (angleDeg(g.normal, v) > 60) continue;
      push(g._rep || ringRep(g, ringChain), ca, 'pc', d); break;
    }
  }
  // ps / ts between rings of opposite chains
  for (const g1 of rings1) {
    for (const g2 of rings2) {
      const v = sub(g2.center, g1.center);
      const d = norm(v);
      const nAng = angleDeg(g1.normal, g2.normal);
      const psi1 = angleDeg(g1.normal, v), psi2 = angleDeg(g2.normal, v);
      if (d <= 7.0 && nAng <= 30 && psi1 <= 45 && psi2 <= 45) {
        push(ringRep(g1, chain1), ringRep(g2, chain2), 'ps', d);
      } else if (d <= 5.0 && Math.abs(90 - nAng) <= 30 && psi1 <= 45 && psi2 <= 45) {
        push(ringRep(g1, chain1), ringRep(g2, chain2), 'ts', d);
      }
    }
  }

  return rows;
}

// representative atom for a ring centroid (for row labels)
function ringRep(group, chain) {
  if (group._rep) return group._rep;
  return { name: 'centroid', resname: 'RNG', resnum: group._resnum || 0, icode: '',
           chain, x: group.center.x, y: group.center.y, z: group.center.z, element: 'C' };
}
