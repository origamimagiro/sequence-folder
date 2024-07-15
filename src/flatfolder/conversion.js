import { M } from "./math.js";
import { CON } from "./constraints.js";
import { NOTE } from "./note.js";
import { AVL } from "./avl.js";
// import { PAR } from "./parallel.js";

export const X = {     // CONVERSION
    L_2_V_EV_EL: (L) => {
        const d = M.min_line_length(L);
        let nV = 0, nE = 0, count = 0;
        const k = 3;    // decrease eps until feature #s repeat sufficiently
        let V_ = [], EV_ = [], EL_ = [], k_ = 1, i_ = 1;
        NOTE.start_check("epsilon");
        for (const i of Array(50).fill().map((_, j) => j + 3)) {
            const eps = d/(2**i);
            NOTE.check(2**i);
            if (eps < M.FLOAT_EPS) { break; }
            try {
                const [V, EV, EL] = X.L_eps_2_V_EV_EL(L, eps);
                if ((V.length == nV) && (EV.length == nE)) {
                    ++count;
                    if (count <= k_) { continue; }
                    V_ = V; EV_ = EV; EL_ = EL; k_ = count; i_ = i;
                    if (count == k) { break; }
                } else {
                    nV = V.length;
                    nE = EV.length;
                    count = 0;
                }
            } catch { nV = 0; nE = 0; count = 0; }
        }
        return [V_, EV_, EL_, i_ - k_];
    },
    L_eps_2_V_EV_EL: (L, eps) => {
        const point_comp = ([x1, y1], [x2, y2]) => {    // point comparator
            if (M.dist([x1, y1], [x2, y2]) < eps) { return 0; }
            const dy = y1 - y2;
            return (Math.abs(dy) >= eps) ? dy : (x1 - x2);
        };
        const line_intersect = ([x1, y1], [x2, y2], [x3, y3], [x4, y4]) => {
            const dx12 = x1 - x2, dx34 = x3 - x4;
            const dy12 = y1 - y2, dy34 = y3 - y4;
            const denom = dx12*dy34 - dx34*dy12;
            if (denom < eps*eps) { return undefined; }
            const x = ((x1*y2 - y1*x2)*dx34 - (x3*y4 - y3*x4)*dx12)/denom;
            const y = ((x1*y2 - y1*x2)*dy34 - (x3*y4 - y3*x4)*dy12)/denom;
            return [x, y];
        };
        const V = [[-Infinity, -Infinity]];     // sentinal first point
        const VL = [[]];                        // lines starting at vertex
        const LV = [];                          // vertex index pair
        const LU = [];                          // unit vector along line
        const LA = [];                          // angle
        const LD = [];                          // distance from origin
        const Q = new AVL((vi, vj) => point_comp(V[vi], V[vj]));
        for (let li = 0; li < L.length; ++li) { // filter endpoints into Q
            let [p, q] = L[li];
            if (point_comp(p, q) > 0) { [p, q] = [q, p]; }
            const [vi, vj] = [p, q].map(v => {
                const vn = V.length;
                V.push(v);
                VL.push([]);
                const j = Q.insert(vn);
                if (j == undefined) { return vn; }
                V.pop();
                VL.pop();
                return j;
            });
            LV.push([vi, vj]);                  // index-based line data
            const u = M.unit(M.sub(V[vj], V[vi]));
            LU.push(u);
            LA.push((u[1] < eps) ? 0 : M.angle(u));
            LD.push(M.dot(M.perp(u), V[vj]));
            VL[vi].push(li);
        }
        const SV = [[undefined, undefined]];    // sentinal segment
        const SU = [[-1, 0]];                   // unit vector along segment
        const SA = [Infinity];                  // angle
        const SD = [undefined];                 // distance from origin
        const SL = [[]];                        // lines overlapping segment
        const point_seg_dist = (vi, si) => SD[si] - M.dot(M.perp(SU[si]), V[vi]);
        const on_line = (vi, si) => (Math.abs(point_seg_dist(vi, si)) < eps);
        let curr;
        const seg_comp = (si, sj) => {
            if(!on_line(curr, si)) {                // assumes si on curr
                throw new Error();
            }
            const dj = point_seg_dist(curr, sj);    // only search near curr,
            if (Math.abs(dj) < eps) {               // so can order locally
                const pi = SV[si][1];
                if ((pi != undefined) && on_line(pi, sj)) { return 0; }
                return SA[sj] - SA[si];
            }
            return -dj;
        };
        const T = new AVL(seg_comp);
        const VP = Array(V).fill();
        const P = [];
        while (Q.length > 0) {
            const vi = Q.remove_next(0);    // processing vertex curr = vi
            curr = vi;
            const v = V[vi];
            const S1 = [];
            SV[0][0] = vi;
            SD[0] = M.dot(M.perp(SU[0]), v);
            if (SA[T.prev(0)] == 0) {
                S1.push(T.remove_prev(0));  // entering horizontal segment
            }
            while (true) {
                const si = T.next(0);
                if ((si == undefined) || !on_line(curr, si)) { break; }
                S1.push(T.remove_next(0));          // other entering segments
            }
            if ((VL[vi].length == 0) && (S1.length < 2)) {
                if (S1.length == 0) { continue; }   // vertex unused
                const si = S1[0];
                let ends = false;
                for (const li of SL[si]) {
                    if (point_comp(V[LV[li][1]], V[vi]) <= 0) {
                        ends = true;    // vertex used (ends some segment)
                        continue;
                    }
                }
                if (!ends) {            // vertex unused (one segment passing)
                    T.insert(si);
                    continue;
                }
            }
            VP[vi] = P.length;
            P.push(v);
            for (const si of S1) {      // close segments
                SV[si][1] = vi;
                for (const li of SL[si]) {
                    if (point_comp(V[LV[li][1]], V[vi]) <= 0) { continue; }
                    VL[vi].push(li);    // add lines that pass through
                }
            }
            VL[vi].sort((i, j) => {     // process lines by distance from vi
                const dj = M.distsq(v, V[LV[j][1]]);
                const di = M.distsq(v, V[LV[i][1]]);
                return dj - di;
            });
            for (const li of VL[vi]) {      // add lines to new segments
                const si = SV.length;
                SV.push([vi, LV[li][1]]);   // endpoint for comparisons
                SU.push(LU[li]);
                SA.push(LA[li]);
                SD.push(LD[li]);
                SL.push([li]);
                const sj = T.insert(si);
                if (sj != undefined) {      // existing segment
                    SV.pop(); SU.pop();
                    SA.pop(); SD.pop();
                    SL.pop();
                    SL[sj].push(li);
                } else {                    // new segment
                    SV[si][1] = undefined;  // erase endpoint
                }
            }
            const pairs = [];
            pairs.push([T.prev(0), T.next(0)]); SA[0] = -Infinity
            pairs.push([T.prev(0), T.next(0)]); SA[0] =  Infinity
            for (const [l, r] of pairs) {   // check adjacent pairs
                if ((l == undefined) || (r == undefined)) {
                    continue;               // incomplete pair
                }
                const vl = SV[l][0], al = SA[l];
                const vr = SV[r][0], ar = SA[r];
                const x = line_intersect(V[vl], M.add(V[vl], SU[l]),
                                         V[vr], M.add(V[vr], SU[r]));
                if ((x == undefined) ||                         // none
                    (point_comp(x, v) == 0) ||                  // near curr
                    ((point_comp(x, v) < 0) && (x[1] <= v[1]))  // behind sweep
                ) { continue; }
                const vx = V.length;        // add point
                V.push(x);
                VL.push([]);
                const vj = Q.insert(vx);
                if (vj != undefined) {      // point near existing point
                    V.pop();
                    VL.pop();
                }
            }
        }
        const X = [];
        const X_map = new Map();
        for (let si = 1; si < SV.length; ++si) {    // combine redundant
            const pp = SV[si].map(vi => VP[vi]);
            if ((pp[0] == undefined) || (pp[1] == undefined)) {
                throw new Error();                  // missing endpoint
            }
            if (pp[1] < pp[0]) { pp.reverse(); }
            const k = M.encode(pp);
            let xi = X_map.get(k);
            if (xi == undefined) {
                xi = X.length;
                X.push([pp, []]);
                X_map.set(k, xi);
            }
            X[xi][1].push(...SL[si]);
        }
        X.sort(([[i1, j1], L1], [[i2, j2], L2]) => {    // lexicographic
            return (i1 == i2) ? (j1 - j2) : (i1 - i2);
        })
        const XP = X.map(xi => xi[0]);
        const XL = X.map(xi => xi[1].sort((a, b) => a - b));
        return [P, XP, XL];
    },
    V_EV_2_VV_FV: (V, EV) => {
        const adj = V.map(() => []);
        for (const [pi, qi] of EV) {
            adj[pi].push(qi);
            adj[qi].push(pi);
        }
        const VV = [];
        for (const [i, v0] of V.entries()) {
            const A = [];
            for (const vi of adj[i]) {
                A.push([vi, M.angle(M.sub(V[vi], v0))]);
            }
            A.sort(([v1, a1], [v2, a2]) => a1 - a2);
            VV.push(A.map(([v, a]) => v));
        }
        const FV = [];
        const seen = new Set();
        for (const [v1, A] of VV.entries()) {
            for (const v2 of A) {
                const key = M.encode([v1, v2]);
                if (!(seen.has(key))) {
                    seen.add(key);
                    const F = [v1];
                    let [i, j] = [v1, v2];
                    while (j != v1) {
                        F.push(j);
                        [i, j] = [j, M.previous_in_list(VV[j], i)];
                        seen.add(M.encode([i, j]));
                    }
                    if (F.length > 2) {
                        FV.push(F);
                    }
                }
            }
        }
        M.sort_faces(FV, V);
        FV.pop(); // remove outer face
        return [VV, FV];
    },
    V_FV_2_VV: (V, FV) => {
        const next = V.map(() => new Map());
        const prev = V.map(() => new Map());
        for (const [fi, V] of FV.entries()) {
            let [v1, v2] = [V[V.length - 2], V[V.length - 1]];
            for (const v3 of V) {
                next[v2].set(v1, v3);
                prev[v2].set(v3, v1);
                [v1, v2] = [v2, v3];
            }
        }
        const VV = V.map(() => []);
        for (const [i, Adj] of VV.entries()) {
            const v0 = next[i].keys().next().value;
            let [v, v_] = [v0, undefined];
            const v1 = prev[i].get(v0);
            if (v1 != undefined) {
                [v, v_] = [v1, prev[i].get(v)];
                while ((v_ != undefined) && (v_ != v0)) {
                    [v, v_] = [v_, prev[i].get(v_)];
                }
            }
            const start = v;
            Adj.push(start);
            v = next[i].get(v);
            while ((v != undefined) && (v != start)) {
                Adj.push(v);
                v = next[i].get(v);
            }
        }
        return VV;
    },
    V_VV_EV_EA_2_VK: (V, VV, EV, EA) => {
        const VVA_map = new Map();
        for (const [i, [v1, v2]] of EV.entries()) {
            const a = EA[i];
            VVA_map.set(M.encode([v1, v2]), a);
            VVA_map.set(M.encode([v2, v1]), a);
        }
        const VK = [];
        for (const [i, A] of VV.entries()) {
            const adj = [];
            let boundary = false;
            let [count_M, count_V, count_U] = [0, 0, 0];
            for (const j of A) {
                const a = VVA_map.get(M.encode([i, j]));
                if (a == "B") {
                    boundary = true;
                    break;
                }
                if (a == "V" || a == "M" || a == "U") {
                    adj.push(j);
                }
                if (a == "M") { ++count_M; }
                if (a == "V") { ++count_V; }
                if (a == "U") { ++count_U; }
            }
            if (boundary || (adj.length == 0)) {
                VK.push(0);
            } else if (
                ((count_U == 0) && (Math.abs(count_M - count_V) != 2)) ||
                (adj.length % 2 != 0)
            ) {                       // violates Maekawa
                VK.push(1);           // far from zero
            } else {
                const angles = adj.map(j => M.angle(M.sub(V[j], V[i])));
                angles.sort((a, b) => a - b);
                let kawasaki = 0;
                for (let j = 0; j < angles.length; j += 2) {
                    kawasaki += angles[j + 1] - angles[j];
                }
                VK.push(Math.abs(kawasaki - Math.PI));
            }
        }
        return VK;
    },
    V_FV_EV_EA_2_Vf_Ff: (V, FV, EV, EA) => {
        const EA_map = new Map();
        for (const [i, vs] of EV.entries()) {
            EA_map.set(M.encode_order_pair(vs), EA[i]);
        }
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const Vf = V.map((p) => undefined);
        const seen = new Set();
        seen.add(0);                    // start search at face 0
        const [v1, v2,] = FV[0];        // first edge of first face
        for (const i of [v1, v2]) {
            Vf[i] = V[i];
        }            // [face, edge, len, parity]
        const Ff = new Array(FV.length);
        const queue = [[0, v1, v2, Infinity, true]];
        let next = 0;
        while (next < queue.length) {                   // Prim's algorithm to
            const [fi, i1, i2, l, s] = queue[next];     // traverse face graph
            Ff[fi] = !s;                                // over spanning tree
            next += 1;                                  // crossing edges of
            const F = FV[fi];                           // maximum length
            const x = M.unit(M.sub(V[i2], V[i1]));
            const y = M.perp(x);
            const xf = M.unit(M.sub(Vf[i2], Vf[i1]));
            const yf = M.perp(xf);
            let vi = F[F.length - 1];
            for (const vj of F) {
                if (Vf[vj] == undefined) {
                    const v = M.sub(V[vj], V[i1]);
                    const dx = M.mul(xf, M.dot(v, x));
                    const dy = M.mul(yf, M.dot(v, y) * (s ? 1 : -1));
                    const vnew = M.add(M.add(dx, dy), Vf[i1]);
                    Vf[vj] = vnew;
                }
                const len = M.distsq(V[vi], V[vj]);
                const f = EF_map.get(M.encode([vi, vj]));
                const a = EA_map.get(M.encode_order_pair([vi, vj]));
                const new_s = (a == "M" || a == "V" || a == "U") ? !s : s;
                if ((f != undefined) && !seen.has(f)) {
                    queue.push([f, vj, vi, len, new_s]);
                    seen.add(f);
                    let prev = len;
                    for (let i = queue.length - 1; i > next; --i) {
                        const curr = queue[i - 1][3];   // O(n^2) but could be
                        if (curr < prev) {              // O(n log n)
                            [queue[i], queue[i - 1]] = [queue[i - 1], queue[i]];
                        } else {
                            break;
                        }
                        prev = curr;
                    }
                }
                vi = vj;
            }
        }
        for (const p of Vf) { if (p == undefined) { debugger; } }
        return [Vf, Ff];
    },
    EV_FV_2_EF_FE: (EV, FV) => {
        const EV_map = new Map();
        for (const [i, V] of EV.entries()) {
            EV_map.set(M.encode(V), i);
        }
        const EF = EV.map(() => [undefined, undefined]);
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                const ei = EV_map.get(M.encode_order_pair([v1, v2]));
                const c = (v2 < v1) ? 0 : 1;
                EF[ei][c] = i;
            }
        }
        for (const [i, F] of EF.entries()) {
            const c = (F[0] == undefined) ? 1 :
                     ((F[1] == undefined) ? 0 : undefined);
            if (c != undefined) {
                EF[i] = [F[c]];
            }
        }
        const FE = FV.map((V) => {
            const E = [];
            let v1 = V[0];
            for (let i = 1; i < V.length; ++i) {
                const v2 = V[i];
                const k = M.encode_order_pair([v1, v2]);
                E.push(EV_map.get(k));
                v1 = v2;
            }
            const k = M.encode_order_pair([v1, V[0]]);
            E.push(EV_map.get(k));
            return E;
        });
        return [EF, FE];
    },
    f_FC_CF_2_fB_set: (f, FC, CF) => {
        const fB = new Set();
        for (const c of FC[f]) {
            for (const f of CF[c]) { fB.add(f); }
        }
        return fB;
    },
    EF_FV_SP_SE_CP_SC_2_CF_FC: (EF, FV, SP, SE, CP, SC) => {
        const SF_map = new Map();
        for (const [i, vs] of SP.entries()) {
            const Fs = [];
            for (const ei of SE[i]) {
                for (const f of EF[ei]) {
                    Fs.push(f);
                }
            }
            SF_map.set(M.encode_order_pair(vs), Fs);
        }
        const SC_map = new Map();
        for (const [i, C] of CP.entries()) {
            let v1 = C[C.length - 1];
            for (const v2 of C) {
                SC_map.set(M.encode([v2, v1]), i);
                v1 = v2;
            }
        }
        const CF = CP.map(() => []);
        const seen = new Set();
        const queue = [];
        for (const [i, Cs] of SC.entries()) {   // Look for a segment on the
            if (Cs.length == 1) {               // border of the overlap graph
                const ci = Cs[0];               // to start BFS
                CF[ci] = SF_map.get(M.encode_order_pair(SP[i]));
                queue.push(ci);
                seen.add(ci);
                break;
            }
        }
        let next = 0;
        while (next < queue.length) {   // BFS on cells in the overlap graph
            const ci = queue[next];
            next++;
            const C = CP[ci];
            let v1 = C[C.length - 1];
            for (const v2 of C) {
                const c = SC_map.get(M.encode([v1, v2]));
                if ((c != undefined) && !seen.has(c)) {
                    queue.push(c);
                    seen.add(c);
                    const Fs = new Set(CF[ci]);
                    const k = M.encode_order_pair([v1, v2]);
                    for (const f of SF_map.get(k)) {
                        const removed = Fs.delete(f);
                        if (!removed) {
                            Fs.add(f);
                        }
                    }
                    CF[c] = Array.from(Fs);
                }
                v1 = v2;
            }
        }
        const FC = FV.map(() => []);
        for (let ci = 0; ci < CF.length; ci++) {
            CF[ci].sort((a, b) => a - b);
            for (const f of CF[ci]) {
                FC[f].push(ci);
            }
        }
        return [CF, FC];
    },
    SE_2_ExE: (SE) => {
        const ExE = new Set();
        for (const edges of SE) {
            for (const [j, v1] of edges.entries()) {
                for (let k = j + 1; k < edges.length; ++k) {
                    const v2 = edges[k];
                    ExE.add(M.encode_order_pair([v1, v2]));
                }
            }
        }
        return Array.from(ExE);
    },
    SE_CF_SC_2_ExF: (SE, CF, SC) => {
        const ExF = new Set();
        for (const [i, C] of SC.entries()) {
            if (C.length == 2) {
                const E = SE[i];
                const [c1, c2] = C;
                const F = [];
                const F1 = new Set(CF[c1]);
                for (const fi of CF[c2]) {
                    if (F1.has(fi)) {
                        F.push(fi);
                    }
                }
                for (const ei of E) {
                    for (const fi of F) {
                        ExF.add(M.encode([ei, fi]));
                    }
                }
            }
        }
        return Array.from(ExF);
    },
    CF_2_BF: (CF) => {                          // O(|C|t^2) <= O(|F|^4)
        const BF_set = new Set();               // t is max faces in a cell
        NOTE.start_check("cell", CF);           // t^2 = O(|B|) <= O(|F|^2)
        for (const [i, F] of CF.entries()) {    // |C| = O(|F|^2)
            NOTE.check(i);
            for (const [j, f1] of F.entries()) {
                for (let k = j + 1; k < F.length; k++) {
                    const f2 = F[k];
                    BF_set.add(M.encode([f1, f2]));
                }
            }
        }
        return Array.from(BF_set).sort();       // |BF| = O(|F|^2)
    },
    EF_SP_SE_CP_CF_2_BF: (EF, SP, SE, CP, CF) => {
        const SF_map = new Map();
        for (const [i, vs] of SP.entries()) {
            const Fs = [];
            for (const ei of SE[i]) {
                for (const f of EF[ei]) {
                    Fs.push(f);
                }
            }
            SF_map.set(M.encode_order_pair(vs), Fs);
        }
        const SC_map = new Map();
        for (const [i, C] of CP.entries()) {
            let v1 = C[C.length - 1];
            for (const v2 of C) {
                SC_map.set(M.encode([v2, v1]), i);
                v1 = v2;
            }
        }
        const BF_set = new Set();
        const seen = new Set();
        const queue = [0];
        {
            const F = CF[0];
            for (let j = 1; j < F.length; ++j) {
                for (let i = 0; i < j; ++i) {
                    BF_set.add(M.encode_order_pair([F[i], F[j]]));
                }
            }
        }
        let next = 0;
        const CF_set = CF.map(F => new Set(F));
        while (next < queue.length) {   // BFS on cells in the overlap graph
            const ci = queue[next++];
            const C = CP[ci];
            let v1 = C[C.length - 1];
            for (const v2 of C) {
                const cj = SC_map.get(M.encode([v1, v2]));
                if ((cj != undefined) && !seen.has(cj)) {
                    queue.push(cj);
                    seen.add(cj);
                    const Fi_set = CF_set[ci];
                    const Fj_set = CF_set[cj];
                    const k = M.encode_order_pair([v1, v2]);
                    for (const fi of SF_map.get(k)) {
                        if (Fi_set.has(fi) || !Fj_set.has(fi)) { continue; }
                        for (const fj of Fj_set) {
                            if (fi == fj) { continue; }
                            BF_set.add(M.encode_order_pair([fi, fj]));
                        }
                    }
                }
                v1 = v2;
            }
        }
        return Array.from(BF_set).sort();
    },
    check_overlap: (p, BI) => {
        return (BI.has(M.encode_order_pair(p)) ? 1 : 0);
    },
    add_constraint: (T, BI, BT) => {
        if (T != undefined) {
            const [type, F] = T;
            const pairs = CON.type_F_2_pairs(type, F);
            for (const p of pairs) {
                const i = BI.get(M.encode_order_pair(p));
                if (i == undefined) { debugger; }
                BT[type][i].push(F);
            }
        }
    },
    BF_BI_EF_ExE_ExF_2_BT0_BT1_BT2: (BF, BI, EF, ExE, ExF) => {
        const BT0 = BF.map(() => []); // taco-taco
        const BT1 = BF.map(() => []); // taco-tortilla
        const BT2 = BF.map(() => []); // tortilla-tortilla
        const BT = [BT0, BT1, BT2];
        NOTE.time("Computing from edge-edge intersections");
        NOTE.start_check("edge-edge intersection", ExE);
        for (const [i, k] of ExE.entries()) {
            NOTE.check(i);
            const [e1, e2] = M.decode(k);
            if ((EF[e1].length != 2) || (EF[e2].length != 2)) { continue; }
            const [f1, f2] = EF[e1];
            const [f3, f4] = EF[e2];
            // note that all of f1, f2, f3, f4 must be unique
            const f1f2 = X.check_overlap([f1, f2], BI);
            const f1f3 = X.check_overlap([f1, f3], BI);
            const f1f4 = X.check_overlap([f1, f4], BI);
            let cons;
            const choice = (f1f2 << 2) | (f1f3 << 1) | f1f4;
            switch (choice) {
                case 0: // 000
                    cons = [CON.T.taco_tortilla, [f3, f4, f2]]; break;
                case 1: // 001
                    cons = [CON.T.tortilla_tortilla, [f1, f2, f4, f3]]; break;
                case 2: // 010
                    cons = [CON.T.tortilla_tortilla, [f1, f2, f3, f4]]; break;
                case 3: // 011
                    cons = [CON.T.taco_tortilla, [f3, f4, f1]]; break;
                case 4: // 100  no overlap
                    break;
                case 5: // 101
                    cons = [CON.T.taco_tortilla, [f1, f2, f4]]; break;
                case 6: // 110
                    cons = [CON.T.taco_tortilla, [f1, f2, f3]]; break;
                case 7: // 111
                    cons = [CON.T.taco_taco, [f1, f2, f3, f4]]; break;
                default: break;
            }
            X.add_constraint(cons, BI, BT);
        }
        NOTE.time("Computing from edge-face intersections");
        NOTE.start_check("edge-face intersection", ExF);
        for (const [i, k] of ExF.entries()) {
            NOTE.check(i);
            const [e, f3] = M.decode(k);
            if (EF[e].length != 2) { continue; }
            const [f1, f2] = EF[e];
            if ((f1 == f3) || (f2 == f3)) { continue; }
            const f1f2 = X.check_overlap([f1, f2], BI);
            let cons;
            if (f1f2 == 1) {
                cons = [CON.T.taco_tortilla, [f1, f2, f3]];
            } else {
                cons = [CON.T.tortilla_tortilla, [f1, f2, f3, f3]];
            }
            X.add_constraint(cons, BI, BT);
        }
        return BT;
    },
    FC_BF_BI_BT0_BT1_2_BT3x: (FC, BF, BI, BT0, BT1) => {
        const BT3x = BF.map(() => new Set());
        NOTE.start_check("variable", BF);
        for (const [i, k] of BF.entries()) {
            NOTE.check(i);
            const [f1, f2] = M.decode(k);
            for (const T of [BT0[i], BT1[i]]) {
                for (const F of T) {
                    for (const f of F) {
                        if ((f == f1) || (f == f2)) { continue; }
                        BT3x[i].add(f);
                    }
                }
            }
        }
        const FG = FC.map(() => new Map());
        NOTE.start_check("taco-tortilla", BT1);
        for (const [i, T] of BT1.entries()) {    // construct connectivity graphs
            NOTE.check(i);
            for (const [a, b, c] of T) {
                const G = FG[c];
                if (!G.has(a)) { G.set(a, []); }
                if (!G.has(b)) { G.set(b, []); }
                G.get(a).push(b);
                G.get(b).push(a);
            }
        }
        const CC = [];
        const seen = new Set();
        NOTE.start_check("face", FG);
        for (const [c, G] of FG.entries()) {    // find connected components
            NOTE.check(c);
            seen.clear();
            for (const [F, _] of G) {
                if (seen.has(F)) { continue; }
                seen.add(F);
                const C = [F];
                let i = 0;
                while (i < C.length) {
                    for (const f of G.get(C[i])) {
                        if (seen.has(f)) { continue; }
                        seen.add(f);
                        C.push(f);
                    }
                    ++i;
                }
                CC.push([c, C]);
            }
            G.clear();
        }
        FG.length = 0;
        NOTE.start_check("component", CC);      // find implied transitivity
        for (const [ci, [c, C]] of CC.entries()) {
            NOTE.check(ci);
            const K = C.map(a => BI.get(M.encode_order_pair([a, c])));
            for (let i = 0; i < C.length - 1; ++i) {
                const a = C[i];
                const kca = K[i];
                if (kca == undefined) { continue; }
                for (let j = i + 1; j < C.length; ++j) {
                    const b = C[j];
                    const kbc = K[j];
                    if (kbc == undefined) { continue; }
                    const kab = BI.get(M.encode_order_pair([a, b]));
                    if (kab == undefined) { continue; }
                    BT3x[kca].add(b);
                    BT3x[kbc].add(a);
                    BT3x[kab].add(c);   // possible a,b,c do not all overlap
                }
            }
            C.length = 0;
            K.length = 0;
            CC[ci] = undefined;
        }
        CC.length = 0;
        for (let i = 0; i < BT3x.length; ++i) {
            BT3x[i] = M.encode(BT3x[i]);
        }
        return BT3x;
    },
    FC_BF_BI_BT0_BT1_W_2_BT3x: async (FC, BF, BI, BT0, BT1, W) => {
        const BT3x = BF.map(() => new Set());
        NOTE.start_check("variable", BF);
        for (const [i, k] of BF.entries()) {
            NOTE.check(i);
            const [f1, f2] = M.decode(k);
            for (const T of [BT0[i], BT1[i]]) {
                for (const F of T) {
                    for (const f of F) {
                        if ((f == f1) || (f == f2)) { continue; }
                        BT3x[i].add(f);
                    }
                }
            }
        }
        const FG = FC.map(() => new Map());
        NOTE.start_check("taco-tortilla", BT1);
        for (const [i, T] of BT1.entries()) {    // construct connectivity graphs
            NOTE.check(i);
            for (const [a, b, c] of T) {
                const G = FG[c];
                if (!G.has(a)) { G.set(a, []); }
                if (!G.has(b)) { G.set(b, []); }
                G.get(a).push(b);
                G.get(b).push(a);
            }
        }
        const CC = [];
        const seen = new Set();
        NOTE.start_check("face", FG);
        for (const [c, G] of FG.entries()) {    // find connected components
            NOTE.check(c);
            seen.clear();
            for (const [F, _] of G) {
                if (seen.has(F)) { continue; }
                seen.add(F);
                const C = [F];
                let i = 0;
                while (i < C.length) {
                    for (const f of G.get(C[i])) {
                        if (seen.has(f)) { continue; }
                        seen.add(f);
                        C.push(f);
                    }
                    ++i;
                }
                CC.push([c, C]);
            }
            G.clear();
        }
        FG.length = 0;
        const P = W.map((_, wi) => new Promise(res => res([wi, undefined])));
        await Promise.all(W.map(w => PAR.send_message(w, "BT3x_start", [BI])));
        const load = 10;
        let i = 0, j = 0;
        NOTE.start_check("component", CC);
        while (j < CC.length) {
            NOTE.check(j);
            const [wi, R] = await Promise.any(P);
            if (R != undefined) {
                for (const [ci, c, cB] of R) {
                    for (const [a, b, kab, kbc, kca] of cB) {
                        BT3x[kbc].add(a);
                        BT3x[kca].add(b);
                        BT3x[kab].add(c);
                    }
                    cB.length = 0;
                    ++j;
                }
                R.length = 0;
            }
            const J = [];
            for (let k = 0; (k < load) && (i < CC.length); ++k, ++i) {
                const [c, C] = CC[i];
                J.push([i, c, C]);
            }
            P[wi] = ((J.length == 0) ? (new Promise(res => {}))
                : PAR.send_message(W[wi], "BT3x", [J]));
            J.length = 0;
        }
        await Promise.all(W.map(w => PAR.send_message(w, "BT3x_stop", [])));
        CC.length = 0;
        for (let i = 0; i < BT3x.length; ++i) {
            BT3x[i] = M.encode(BT3x[i]);
        }
        return BT3x;
    },
    FC_CF_BF_BT3x_2_BT3: (FC, CF, BF, BT3x) => { // O(|B|kt) <= O(|F|^5)
        const BT3 = [];                          // k is max cells in a face,
        const FC_sets = FC.map(C => new Set(C)); // t is max faces in a cell
        let nx = 0;                              // k = O(|C|) <= O(|F|^2)
        NOTE.start_check("variable", BF);        // t = O(|F|)
        for (const [i, k] of BF.entries()) {     // |B| = O(|F|^2)
            NOTE.check(i);
            const [f1, f2] = M.decode(k);
            const C = FC_sets[f1];
            const S = new Set(M.decode(BT3x[i]));
            delete BT3x[i];
            const X = new Set();
            const T = new Set();
            for (const c of FC[f2]) {
                if (C.has(c)) {
                    for (const f3 of CF[c]) {
                        if ((f3 == f1) || (f3 == f2)) { continue; }
                        if (S.has(f3)) { X.add(f3); }
                        else           { T.add(f3); }
                    }
                }
            }
            nx += X.size;
            BT3.push(M.encode(T));
        }
        return [BT3, nx/3];     // |BT3| = O(|B||F|) <= O(|F|^3)
    },
    FC_CF_BF_BT3x_W_2_BT3: async (FC, CF, BF, BT3x, W) => {
        const P = W.map((_, wi) => new Promise(res => res([wi, undefined])));
        await Promise.all(W.map(w => PAR.send_message(w, "BT3_start", [FC, CF])));
        const BT3 = BF.map(() => undefined);
        const load = 100;
        let i = 0, j = 0, nx = 0;
        NOTE.start_check("variable", BF);
        while (j < BF.length) {
            NOTE.check(j);
            const [wi, R] = await Promise.any(P);
            if (R != undefined) {
                for (const [bi, bF, nxi] of R) {
                    BT3[bi] = bF;
                    ++j;
                    nx += nxi;
                }
            }
            const J = [];
            for (let k = 0; (k < load) && (i < BF.length); ++k, ++i) {
                J.push([i, BF[i], BT3x[i]]);
                delete BT3x[i];
            }
            P[wi] = ((J.length == 0) ? (new Promise(res => {}))
                : PAR.send_message(W[wi], "BT3", [J]));
            J.length = 0;
        }
        await Promise.all(W.map(w => PAR.send_message(w, "BT3_stop", [])));
        return [BT3, nx/3];
    },
    EF_SP_SE_CP_FC_CF_BF_BT3x_2_BT3: (EF, SP, SE, CP, FC, CF, BF, BT3x) => {
        const FC_set = FC.map(C => new Set(C));
        const CF_set = CF.map(F => new Set(F));
        const SC_map = new Map();
        for (const [i, C] of CP.entries()) {
            let v1 = C[C.length - 1];
            for (const v2 of C) {
                SC_map.set(M.encode([v2, v1]), i);
                v1 = v2;
            }
        }
        const SF_map = new Map();
        for (const [i, [v1, v2]] of SP.entries()) {
            const F1 = [];
            const F2 = [];
            for (const ei of SE[i]) {
                for (const f of EF[ei]) {
                    const c1 = SC_map.get(M.encode([v1, v2]));
                    const c2 = SC_map.get(M.encode([v2, v1]));
                    if ((c1 != undefined) && CF_set[c1].has(f)) { F1.push(f); }
                    if ((c2 != undefined) && CF_set[c2].has(f)) { F2.push(f); }
                }
            }
            SF_map.set(M.encode([v1, v2]), F1);
            SF_map.set(M.encode([v2, v1]), F2);
        }
        const BT3 = [];
        let nx = 0;
        const C = new Set();
        const X = new Set();
        const T = new Set();
        NOTE.start_check("variable", BF);
        for (const [i, k] of BF.entries()) {
            NOTE.check(i);
            let [f1, f2] = M.decode(k);
            if (FC[f1].length > FC[f2].length) { [f1, f2] = [f2, f1]; }
            C.clear();
            X.clear();
            T.clear();
            for (const ci of FC[f1]) {
                if (FC_set[f2].has(ci)) {
                    C.add(ci);
                }
            }
            const S = new Set(M.decode(BT3x[i]));
            delete BT3x[i];
            const seen = new Set();
            for (const ci of C) {
                if (seen.has(ci)) { continue; }
                for (const f3 of CF[ci]) {
                    if ((f3 == f1) || (f3 == f2)) { continue; }
                    if (S.has(f3)) { X.add(f3); }
                    else           { T.add(f3); }
                }
                const Q = [ci];
                seen.add(ci);
                let qi = 0;
                while (qi < Q.length) {
                    const cj = Q[qi++];
                    const P = CP[cj];
                    let v1 = P[P.length - 1];
                    for (const v2 of P) {
                        const cj = SC_map.get(M.encode([v1, v2]));
                        if ((cj != undefined) && C.has(cj) && !seen.has(cj)) {
                            Q.push(cj);
                            seen.add(cj);
                            const k = M.encode([v1, v2]);
                            for (const f3 of SF_map.get(k)) {
                                if ((f3 == f1) || (f3 == f2)) { continue; }
                                if (S.has(f3)) { X.add(f3); }
                                else           { T.add(f3); }
                            }
                        }
                        v1 = v2;
                    }
                }
            }
            nx += X.size;
            BT3.push(M.encode(T));
        }
        return [BT3, nx/3];
    },
    BF_GB_GA_GI_2_edges: (BF, GB, GA, GI) => {
        const edges = [];
        for (const [i, B] of GB.entries()) {
            const orders = M.bit_decode(GA[i][GI[i]], B.length);
            for (const [j, F] of B.entries()) {
                const [f1, f2] = M.decode(BF[F]);
                const o = orders[j];
                edges.push(M.encode((o == 1) ? [f1, f2] : [f2, f1]));
            }
        }
        return edges;
    },
    edges_Ff_2_FO: (edges, Ff) => {
        return edges.map((k) => {
            const [f1, f2] = M.decode(k);
            return [f1, f2, (Ff[f2] ? 1 : -1)];
        });
    },
    CF_edges_2_CD: (CF, edges) => {
        const edge_map = new Set(edges);
        return CF.map((F, i) => {
            const S = F.map(i => i);
            S.sort((a, b) => (edge_map.has(M.encode([a, b])) ? 1 : -1));
            return S;
        });
    },
    Ctop_SC_SE_EF_Ff_2_SD: (Ctop, SC, SE, EF, Ff) => {
        const EF_set = new Set(
            EF.filter(F => F.length == 2).map(F => M.encode_order_pair(F)));
        const SD = SC.map((C, si) => {
            const F = C.map(ci => Ctop[ci]);
            if (F[0] == F[1]) { return "N"; }
            if ((F[0] == undefined) || (F[1] == undefined)) { return "B"; }
            const flips = F.map(fi => Ff[fi]);
            if ((flips[0] == flips[1]) &&
                EF_set.has(M.encode_order_pair(F))) { return "C"; }
            let left = false, right = false;
            for (const ei of SE[si]) {
                const [fi, fj] = EF[ei];
                if (Ff[fi] == Ff[fj]) { continue; }
                if ((fi == F[0]) || (fj == F[0])) { left = true; }
                if ((fi == F[1]) || (fj == F[1])) { right = true; }
            }
            if (left == right) { return "B"; }
            return left ? "BL" : "BR";
        });
        return SD;
    },
    Ctop_CP_SC_SD_Ff_P_2_RP_Rf: (Ctop, CP, SC, SD, Ff, P) => {
        const cn = Ctop.length;
        const CC = Array(cn).fill(0).map(() => []);
        for (let si = 0; si < SC.length; ++si) {
            const C = SC[si];
            if (C.length != 2) { continue; }
            const d = SD[si];
            if (SD[si][0] == "B") { continue; }
            const [ci, cj] = C;
            CC[ci].push(cj);
            CC[cj].push(ci);
        }
        const RP = [];
        const Rf = [];
        let ri = 0;
        const seen = Array(cn).fill(false);
        for (let ci = 0; ci < cn; ++ci) {
            const fi = Ctop[ci];
            if ((fi == undefined) || seen[ci]) { continue; }
            seen[ci] = true;
            const C = [ci];
            let i = 0;
            while (i < C.length) {
                for (const cj of CC[C[i]]) {
                    if (!seen[cj]) {
                        seen[cj] = true;
                        C.push(cj);
                    }
                }
                ++i;
            }
            const Adj = P.map(() => new Set());
            for (const ci of C) {
                const P_ = CP[ci];
                let pi = P_[P_.length - 1];
                for (const pj of P_) {
                    const A = Adj[pj];
                    if (A.has(pi)) {
                        A.delete(pi);
                    } else {
                        Adj[pi].add(pj);
                    }
                    pi = pj;
                }
            }
            let start;
            for (let i = 0; i < Adj.length; ++i) {
                if (Adj[i].size == 1) { start = i; break; }
            }
            const Q = [];
            let u = start;
            do {
                Q.push(u);
                let v;
                for (v of Adj[u]) { break; }
                u = v;
            } while (u != start);
            RP.push(Q);
            Rf.push(Ff[fi]);
        }
        return [RP, Rf];
    },
};
