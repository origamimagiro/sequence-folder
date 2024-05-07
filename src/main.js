import { M } from "./flatfolder/math.js";
import { NOTE } from "./flatfolder/note.js";
import { SVG } from "./flatfolder/svg.js";
import { IO } from "./flatfolder/io.js";
import { X } from "./flatfolder/conversion.js";
import { SOLVER } from "./flatfolder/solver.js";
import { CON } from "./flatfolder/constraints.js";

window.onload = () => { MAIN.startup(); };  // entry point

const MAIN = {
    color: {
        background: "lightgray",
        face: {
            top: "#AAA",
            bottom: "#FFF",
        },
        edge: {
            U: "black",
            F: "lightgray",
            B: "black",
        },
        rand: [
            "lightpink", "lightgreen", "lightskyblue", "gold",
            "lightsalmon", "powderblue", "lavender", "sandybrown"
        ],
    },
    startup: () => {
        CON.build();
        NOTE.clear_log();
        NOTE.start("*** Starting Flat-Folder ***");
        NOTE.time("Initializing interface");
        const [b, s] = [50, SVG.SCALE];
        const main = document.getElementById("main");
        for (const [k, v] of Object.entries({
            xmlns: SVG.NS,
            style: "background: lightgray",
            viewBox: [0, 0, 3*s, s].join(" "),
        })) {
            main.setAttribute(k, v);
        }
        document.getElementById("shadow").value = 0;
        for (const [i, id] of ["input", "cp", "output"].entries()) {
            const svg = document.getElementById(id);
            for (const [k, v] of Object.entries({
                xmlns: SVG.NS,
                height: s,
                width: s,
                x: i*s,
                y: 0,
                viewBox: [-b, -b, s + 2*b, s + 2*b].join(" "),
            })) {
                svg.setAttribute(k, v);
            }
        }
        document.getElementById("import").onchange = (e) => {
            if (e.target.files.length > 0) {
                const file_reader = new FileReader();
                file_reader.onload = MAIN.process_file;
                file_reader.readAsText(e.target.files[0]);
            }
        };
    },
    process_file: (e) => {
        NOTE.clear_log();
        NOTE.start("*** Starting File Import ***");
        const doc = e.target.result;
        const file_name = document.getElementById("import").value;
        const parts = file_name.split(".");
        const type = parts[parts.length - 1].toLowerCase();
        if (type != "fold") {
            console.log(`Found file with extension ${type}, FOLD format required`);
            return;
        }
        NOTE.time(`Importing from file ${file_name}`);
        const FILE = JSON.parse(doc);
        if ((FILE.file_frames == undefined) || (FILE.file_frames.length < 2)) {
            console.log("File does not have at least 2 FOLD frames to sequence");
            return;
        }
        FILE.i = 0;
        MAIN.draw_frame(FILE);
        document.getElementById("next").onclick = () => {
            FILE.i = Math.min(FILE.i + 1, FILE.file_frames.length - 2);
            MAIN.draw_frame(FILE);
        };
        document.getElementById("prev").onclick = () => {
            FILE.i = Math.max(FILE.i - 1, 0);
            MAIN.draw_frame(FILE);
        };
        document.getElementById("flip").onchange = () => {
            MAIN.draw_frame(FILE);
        };
        document.getElementById("shadow").onchange = () => {
            MAIN.draw_frame(FILE);
        };
        // console.log(FILE);
    },
    get_frame: (FILE, i) => {
        const frame = FILE.file_frames[i];
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(
            frame.vertices_coords,
            frame.faces_vertices
        );
        FOLD.FL = frame["faces_lf:group"];
        FOLD.FO = frame.faceOrders;
        FOLD.line = frame["lf:line"];
        FOLD.points = frame["lf:points"];
        const edges = FOLD.FO.map(([f1, f2, o]) => {
            return M.encode(((FOLD.Ff[f2] ? 1 : -1)*o >= 0) ? [f1, f2] : [f2, f1]);
        });
        CELL.CD = X.CF_edges_2_CD(CELL.CF, edges);
        FOLD.EA = MAIN.EF_Ff_edges_2_EA(FOLD.EF, FOLD.Ff, edges);
        FOLD.Vf = X.V_FV_EV_EA_2_Vf_Ff(FOLD.V, FOLD.FV, FOLD.EV, FOLD.EA)[0];
        if (M.polygon_area2(M.expand(FOLD.FV[0], FOLD.Vf)) < 0) {
            FOLD.Vf = FOLD.Vf.map(v => M.add(M.refY(v), [0, 1]));
        }
        const v0 = FOLD.Vf[0];
        FOLD.Vf = FOLD.Vf.map(p => M.sub(p, v0));
        const [c1, s1] = FOLD.Vf[1];
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, c1, -s1));
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, 0, 1));
        FOLD.Vf = M.normalize_points(FOLD.Vf);
        return [FOLD, CELL];
    },
    draw_frame: (FILE) => {
        const [F1, C1] = MAIN.get_frame(FILE, FILE.i);
        if (FILE.i < FILE.file_frames.length - 1) {
            const out = SVG.clear("output");
            const [F2, C2] = MAIN.get_frame(FILE, FILE.i + 1);
            MAIN.draw_cp(SVG.clear("cp"), F2);
            MAIN.draw_state(out, F2, C2);
            MAIN.draw_state(SVG.clear("input"), F1, C1, F2);
        } else {
            MAIN.draw_state(SVG.clear("input"), F1, C1);
            MAIN.draw_cp(SVG.clear("cp"), F1, false);
        }
    },
    EF_Ff_edges_2_EA: (EF, Ff, edges) => {
        const edge_map = new Set(edges);
        return EF.map(F => {
            if (F.length != 2) { return "B"; }
            const [i, j] = F;
            if (edge_map.has(M.encode([i, j]))) { return Ff[i] ? "M" : "V"; }
            if (edge_map.has(M.encode([j, i]))) { return Ff[i] ? "V" : "M"; }
            return "F";
        });
    },
    draw_cp: (svg, F, bold = true) => {
        // console.log(F);
        const {V, Vf, FV, EV, EF, EA, Ff, FO, FL} = F;
        const faces = FV.map(F => M.expand(F, Vf));
        const lines = EV.map(E => M.expand(E, Vf));
        const colors = EA.map(a => {
            if (a == "B") { return "black"; }
            if (a == "M") { return "blue"; }
            if (a == "V") { return "red"; }
            if (a == "F") { return "gray"; }
        });
        const g1 = SVG.append("g", svg, {id: "flat_f"});
        const g2 = SVG.append("g", svg, {id: "flat_e"});
        const g3 = SVG.append("g", svg, {id: "flat_p"});
        if ((FL == undefined) || !bold) {
            SVG.draw_segments(g2, lines, {stroke: colors, id: true});
            SVG.draw_polygons(g1, faces, {fill: "white", id: true});
        } else {
            const Fc = FL.map(i => (i < 0) ? "white" : (
                MAIN.color.rand[i % MAIN.color.rand.length]));
            SVG.draw_polygons(g1, faces, {fill: Fc, id: true});
            SVG.draw_segments(g2, lines, {stroke: colors, id: true, stroke_width: 1,
                filter: (i) => {
                    const [f, g] = EF[i];
                    return (g == undefined) || (FL[f] == FL[g]);
                },
            });
            SVG.draw_segments(g2, lines, {stroke: colors, id: true, stroke_width: 5,
                filter: (i) => {
                    const [f, g] = EF[i];
                    return (g != undefined) && (FL[f] != FL[g]);
                },
            });
        }
        // SVG.draw_points(g3, Vf, {text: true, fill: "green"});
    },
    draw_state: (svg, FOLD, CELL, F2) => {
        const {Ff, EF, FO} = FOLD;
        const {P, CP, CD} = CELL;
        const flip = document.getElementById("flip").checked;
        const m = [0.5, 0.5];
        const Q = M.normalize_points(
            P.map(p => (flip ? M.add(M.refX(M.sub(p, m)), m) : p))
        );
        const Ctop = CD.map(S => flip ? S[0] : S[S.length - 1]);
        const [UP, UF, SP, SD] = X.tops_CP_EF_Ff_P_2_UP_UF_SP_SD(Ctop, CP, EF, Ff, Q);
        const Ucolor = UF.map(d => {
            if (d == undefined) { return undefined; }
            if (Ff[d] != flip)  { return MAIN.color.face.top; }
            else                { return MAIN.color.face.bottom; }
        });
        const cells = UP.map(V => M.expand(V, Q));
        const G = {};
        for (const id of ["c", "shadow", "s_crease", "s_edge"]) {
            G[id] = SVG.append("g", svg, {id: `${svg.id}_${id}`});
        }
        SVG.draw_polygons(G.c, cells, {
            id: true, fill: Ucolor, stroke: "none"});
        const shadow = +document.getElementById("shadow").value;
        if (shadow > 0) {
            SVG.draw_shadows(G.shadow, cells, EF, Ff, CD, UP, UF, Q, flip, shadow);
        }
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(G.s_crease, lines, {
            id: true, stroke: MAIN.color.edge.F,
            filter: (i) => SD[i] == "C"});
        SVG.draw_segments(G.s_edge, lines, {
            id: true, stroke: MAIN.color.edge.B,
            filter: (i) => SD[i] == "B"});
        if ((F2 != undefined) && (F2.points != undefined)) {
            const line = [MAIN.line_2_coords(F2.line).map(
                p => flip ? M.add(M.refX(M.sub(p, m)), m): p
            )];
            SVG.draw_segments(G.s_edge, line, {
                id: true, stroke: "purple", stroke_width: 5,}
            );
            SVG.draw_points(G.s_edge, M.expand(F2.points, Q),
                {fill: "green", r: 10}
            );
        }
    },
    line_2_coords: (line) => {
        const [u, d] = line;
        const p = M.mul(u, d);
        const off = M.mul(M.perp(u), 10);
        const p1 = M.add(p, off);
        const p2 = M.sub(p, off);
        return [p1, p2];
    },
    V_FV_2_FOLD_CELL: (V, FV) => {
        const Ff = MAIN.FV_V_2_Ff(FV, V);
        const EV_set = new Set();
        for (const fV of FV) {
            let i = fV.length - 1;
            for (let j = 0; j < fV.length; ++j) {
                EV_set.add(M.encode_order_pair([fV[i], fV[j]]));
                i = j;
            }
        }
        const EV = Array.from(EV_set).sort().map(k => M.decode(k));
        const [EF, FE] = X.EV_FV_2_EF_FE(EV, FV);
        const L = EV.map(vs => vs.map(i => V[i]));
        const eps = M.min_line_length(L) / M.EPS;
        NOTE.time(`Using eps ${eps} from min line length ${
            eps*M.EPS} (factor ${M.EPS})`);
        NOTE.time("Constructing points and segments from edges");
        const [P, SP, SE] = X.L_2_V_EV_EL(L, eps);
        NOTE.annotate(P, "points_coords");
        NOTE.annotate(SP, "segments_points");
        NOTE.annotate(SE, "segments_edges");
        NOTE.lap();
        NOTE.time("Constructing cells from segments");
        const [PP,CP] = X.V_EV_2_VV_FV(P, SP);
        NOTE.annotate(CP, "cells_points");
        NOTE.lap();
        NOTE.time("Computing segments_cells");
        const [SC, CS] = X.EV_FV_2_EF_FE(SP, CP);
        NOTE.annotate(SC, "segments_cells");
        NOTE.annotate(CS, "cells_segments");
        NOTE.lap();
        NOTE.time("Making face-cell maps");
        const [CF, FC] = X.EF_FV_SP_SE_CP_SC_2_CF_FC(EF, FV, SP, SE, CP, SC);
        const BF = X.CF_2_BF(CF);
        NOTE.annotate(BF, "variables_faces");
        NOTE.lap();
        const FOLD = {V, FV, EV, EF, FE, Ff, eps};
        const CELL = {P, SP, SE, PP, CP, CS, SC, CF, FC, BF};
        return [FOLD, CELL];
    },
    linearize: (edges, n) => {
        const Adj = Array(n).fill(0).map(() => []);
        for (const s of edges) {
            const [f1, f2] = M.decode(s);
            Adj[f1].push(f2);
        }
        const L = [];
        const seen = Array(n).fill(false);
        const dfs = (i) => {
            if (seen[i]) { return; }
            seen[i] = true;
            for (const j of Adj[i]) {
                dfs(j);
            }
            L.push(i);
        };
        for (let i = 0; i < n; ++i) {
            dfs(i);
        }
        L.reverse();
        console.assert(L.length == n);
        const idx_map = Array(n).fill(undefined);
        for (let i = 0; i < n; ++i) {
            const fi = L[i];
            idx_map[fi] = i;
        }
        for (const s of edges) {
            const [f1, f2] = M.decode(s);
            if (idx_map[f1] > idx_map[f2]) {
                return undefined; // cycle
            }
        }
        for (let i = 0; i < n; ++i) {
            seen[i] = false;
        }
        const layers = [];
        for (let i = 0; i < n; ++i) {
            const fi = L[i];
            if (seen[fi]) { continue; }
            seen[fi] = true;
            const layer = [fi];
            const Adj_set = new Set();
            for (const fj of Adj[fi]) {
                Adj_set.add(fj);
            }
            for (let j = i + 1; j < L.length; ++j) {
                const fj = L[j];
                if (seen[fj]) { continue; }
                if (!Adj_set.has(fj)) {
                    seen[fj] = true;
                    layer.push(fj);
                }
                for (const fk of Adj[fj]) {
                    Adj_set.add(fk);
                }
            }
            layers.push(layer);
        }
        return layers;
    },
    FV_V_2_Ff: (FV, V) => FV.map(fV => (M.polygon_area2(fV.map(i => V[i])) < 0)),
};
