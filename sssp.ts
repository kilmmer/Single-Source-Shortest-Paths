// Implementation of the Single-Source Shortest Paths (SSSP) algorithm described in the paper "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"
// Authors: Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin (arXiv:2504.17033v2, 30 Jul 2025)

// Main decision: The algorithm is implemented in TypeScript to be executable in Node.js environments. We assume directed graphs with non-negative real weights.
// The graph is represented as a Map<number, {to: number, w: number}[]> for adjacencies, vertices from 0 to n-1.
// To simplify, we ignore the transformation to constant degrees, as the time remains correct for general graphs (the article uses it for asymptotic analysis).
// We assume that there are no ties in numerical distances, as per the article's assumption (Assumption 2.1). If there are ties, the behavior may be incorrect, but with real weights, it is rare.
// We use extended comparison for relaxations to simulate total order in paths, but priorities in heaps and structures use only the numerical value (for simplicity and fidelity).
// Parameters k and t are calculated with base 2 logarithm, as common in complexity.
// Infinity is used for initial distances; care with precision in floating point (not treated, as comparison-addition model).
// The implementation is recursive for BMSSP, which may limit depth in large n (depth O(log^{1/3} n) ~ small).
// Reflection: The idea of the article is innovative, breaking the sorting barrier with boundary reduction, but in practice, for small n, Dijkstra with heap is simpler and faster. For large n, this implementation can be tested, but the constant overhead is high due to the complex data structure.
// Counterpoint: If the graph has negative weights, the algorithm fails (as stated in the article). Also, the deterministic nature is good, but randomized algorithms like in [DMSY23] may be faster on average.
// Alternative perspective: For sparse graphs, Thorup's algorithm for integer weights is linear, but here it is for reals, so it is useful in theoretical cases.

// Step 1: Definition of auxiliary structures

interface Graph {
    n: number;
    adj: Map<number, { to: number, w: number }[]>;
}

// Class for Priority Queue with decreaseKey for base case (binary heap)
class PriorityQueue {
    private heap: { v: number, priority: number }[] = [];
    private position: Map<number, number> = new Map();

    public has(v: number): boolean {
        return this.position.has(v);
    }

    insert(v: number, priority: number) {
        const idx = this.heap.length;
        this.heap.push({ v, priority });
        this.position.set(v, idx);
        this.siftUp(idx);
    }

    extractMin(): number | null {
        if (this.heap.length === 0) return null;
        const first = this.heap[0];
        if (!first) return null;
        const min = first.v;
        const last = this.heap.pop()!;
        if (this.heap.length > 0) {
            this.heap[0] = last;
            this.position.set(last.v, 0);
            this.siftDown(0);
        }
        this.position.delete(min);
        return min;
    }

    decreaseKey(v: number, priority: number) {
        const idx = this.position.get(v);
        if (idx === undefined || this.heap[idx] == null || priority > this.heap[idx].priority) return;
        this.heap[idx].priority = priority;
        this.siftUp(idx);
    }

    isEmpty(): boolean {
        return this.heap.length === 0;
    }

    private siftUp(i: number) {
        while (i > 0) {
            const p = Math.floor((i - 1) / 2);
            if (
                this.heap[i] != null &&
                this.heap[p] != null &&
                this.heap[i]!.priority < this.heap[p]!.priority
            ) {
                this.swap(i, p);
                i = p;
            } else break;
        }
    }

    private siftDown(i: number) {
        const size = this.heap.length;
        while (true) {
            let c = 2 * i + 1;
            if (c >= size) break;
            if (
                c + 1 < size &&
                this.heap[c + 1] != null &&
                this.heap[c] != null &&
                this.heap[c + 1]!.priority < this.heap[c]!.priority
            ) c++;
            if (this.heap[c]!.priority < this.heap[i]!.priority) {
                this.swap(i, c);
                i = c;
            } else break;
        }
    }

    private swap(i: number, j: number) {
        const temp = this.heap[i]!;
        this.heap[i] = this.heap[j]!;
        this.heap[j] = temp;
        this.position.set(this.heap[i].v, i);
        this.position.set(this.heap[j].v, j);
    }
}

// Data structure for Lemma 3.3 (block-based for partial sorting)
// Decision: Blocks as arrays for simplicity. Location map for fast deletes. Binary search in D1 as number of blocks is small.
// Split uses sort O(M log M), amortized fine. Pull uses sort O( collected log collected ) with collected O(M).
class PartialSortDS {
    private M: number;
    private B: number;
    private D0: { key: number, value: number }[][] = [];
    private D1: { block: { key: number, value: number }[], upper: number }[] = [];
    private location: Map<number, { d: 0 | 1, blockIdx: number, itemIdx: number }> = new Map();

    constructor(M: number, B: number) {
        this.M = M;
        this.B = B;
        this.D1.push({ block: [], upper: B });
    }

    private _delete(key: number) {
        const loc = this.location.get(key);
        if (!loc) return;
        if (loc.d === 0) {
            this.D0[loc.blockIdx]!.splice(loc.itemIdx, 1);
            // Update itemIdx for remaining in block
            for (let j = loc.itemIdx; j < this.D0[loc.blockIdx]!.length; j++) {
                this.location.get(this.D0[loc.blockIdx]![j]!.key)!.itemIdx = j;
            }
            if (this.D0[loc.blockIdx]!.length === 0) {
                this.D0.splice(loc.blockIdx, 1);
                // Update blockIdx for subsequent blocks
                for (let b = loc.blockIdx; b < this.D0.length; b++) {
                    this.D0[b]!.forEach((item, j) => this.location.get(item.key)!.blockIdx = b);
                }
            }
        } else {
            this.D1[loc.blockIdx]!.block.splice(loc.itemIdx, 1);
            for (let j = loc.itemIdx; j < this.D1[loc.blockIdx]!.block.length; j++) {
                this.location.get(this.D1[loc.blockIdx]!.block[j]!.key)!.itemIdx = j;
            }
            if (this.D1[loc.blockIdx]!.block.length === 0) {
                this.D1.splice(loc.blockIdx, 1);
                for (let b = loc.blockIdx; b < this.D1.length; b++) {
                    this.D1[b]!.block.forEach((item, j) => this.location.get(item.key)!.blockIdx = b);
                }
            }
        }
        this.location.delete(key);
    }

    insert(key: number, value: number) {
        const currentValue = this.location.has(key) ? this.getValue(key) : Infinity;
        if (value >= currentValue) return;
        this._delete(key);
        // Binary search for block in D1
        let left = 0, right = this.D1.length - 1;
        while (left <= right) {
            const mid = Math.floor((left + right) / 2);
            if (this.D1[mid]!.upper >= value) right = mid - 1;
            else left = mid + 1;
        }
        let blockIdx = left;
        if (blockIdx > this.D1.length - 1) blockIdx = this.D1.length - 1;
        const block = this.D1[blockIdx]!.block;
        block.push({ key, value });
        this.location.set(key, { d: 1, blockIdx, itemIdx: block.length - 1 });
        if (block.length > this.M) this._splitD1(blockIdx);
    }

    private _splitD1(blockIdx: number) {
        let block = this.D1[blockIdx]!.block;
        block.sort((a, b) => a.value - b.value);
        const medianIdx = Math.floor(block.length / 2);
        const median = block[medianIdx]!.value;
        const block1 = block.filter(item => item.value < median);
        const block2 = block.filter(item => item.value >= median);
        let upper1 = block1.length > 0 ? block1.reduce((max, item) => Math.max(max, item.value), -Infinity) : this.D1[blockIdx]!.upper;
        let upper2 = block2.length > 0 ? block2.reduce((max, item) => Math.max(max, item.value), -Infinity) : this.D1[blockIdx]!.upper;
        this.D1.splice(blockIdx, 1, { block: block1, upper: upper1 }, { block: block2, upper: upper2 });
        block1.forEach((item, j) => this.location.set(item.key, { d: 1, blockIdx, itemIdx: j }));
        block2.forEach((item, j) => this.location.set(item.key, { d: 1, blockIdx: blockIdx + 1, itemIdx: j }));
        for (let b = blockIdx + 2; b < this.D1.length; b++) {
            this.D1[b]!.block.forEach((item, j) => this.location.set(item.key, { d: 1, blockIdx: b, itemIdx: j }));
        }
    }

    batchPrepend(L: { key: number, value: number }[]) {
        const map: Map<number, number> = new Map<number, number>();
        L.forEach(({ key, value }) => map.set(key, Math.min(map.get(key) ?? Infinity, value)));
        const newL: Array<any> = [];
        for (let [key, value] of map) {
            const current = this.location.has(key) ? this.getValue(key) : Infinity;
            if (value >= current) continue;
            this._delete(key);
            newL.push({ key, value });
        }
        if (newL.length === 0) return;
        newL.sort((a, b) => a.value - b.value);
        const blockSize = Math.ceil(this.M / 2);
        const newBlocks: Array<any> = [];
        for (let i = 0; i < newL.length; i += blockSize) {
            newBlocks.push(newL.slice(i, i + blockSize));
        }
        this.D0 = newBlocks.concat(this.D0);
        this.D0.forEach((block, b) => block.forEach((item, j) => this.location.set(item.key, { d: 0, blockIdx: b, itemIdx: j })));
    }

    pull() {
        let collected: { key: number, value: number }[] = [];
        let d0Count = 0;
        for (let b = 0; b < this.D0.length; b++) {
            collected.push(...this.D0[b]!);
            d0Count += this.D0[b]!.length;
            if (collected.length > this.M) break;
        }
        let d1Start = 0;
        for (let b = 0; b < this.D1.length; b++) {
            collected.push(...this.D1[b]!.block);
            d1Start = b + 1;
            if (collected.length > this.M) break;
        }
        if (collected.length <= this.M) {
            const S = collected.map(item => item.key);
            this.D0 = [];
            this.D1 = [{ block: [], upper: this.B }];
            this.location.clear();
            return { x: this.B, S };
        }
        collected.sort((a, b) => a.value - b.value);
        const S_items = collected.slice(0, this.M);
        const S = S_items.map(item => item.key);
        const x = collected[this.M]!.value; // Approximation for min remaining, assuming unique
        // Delete S
        const SSet = new Set(S);
        this.D0 = this.D0.filter((block, b) => {
            block = block.filter(item => !SSet.has(item.key));
            if (block.length === 0) return false;
            block.forEach((item, j) => this.location.set(item.key, { d: 0, blockIdx: b, itemIdx: j }));
            return true;
        });
        this.D1 = this.D1.filter((item, b) => {
            item.block = item.block.filter(it => !SSet.has(it.key));
            if (item.block.length === 0) return false;
            item.block.forEach((it, j) => this.location.set(it.key, { d: 1, blockIdx: b, itemIdx: j }));
            return true;
        });
        return { x, S };
    }

    getValue(key: number): number {
        const loc = this.location.get(key);
        if (!loc) return Infinity;
        if (loc.d === 0) return this.D0[loc.blockIdx]![loc.itemIdx]!.value;
        return this.D1[loc.blockIdx]!.block[loc.itemIdx]!.value;
    }

    isEmpty() {
        return this.location.size === 0;
    }

}

// Global variable to enable/disable logs
const enableLogs = !process.argv.includes('--no-log');
const logHistory: string[] = []; // Store log history
let logCounter = 0; // Global counter for log entries
const maxTotalLogs = 1000; // Maximum number of logs allowed during execution

function log(...args: any[]) {
    if (logCounter >= maxTotalLogs) return; // Stop logging if the limit is reached
    const logMessage = `[Log ${logCounter++}] ` + args.join(' ');
    logHistory.push(logMessage);
    if (logHistory.length > 10) {
        logHistory.shift(); // Keep only the last 10 logs
    }
    if (enableLogs) {
        console.log(logMessage);
    }
}

function displayLogHistory() {
    console.log('Last 10 logs:');
    logHistory.forEach((message) => console.log(message));
}

// Step 2: Auxiliary functions for the algorithm
function findPivots(B: number, S: number[], d: number[], depth: number[], Pred: number[], adj: Map<number, { to: number, w: number }[]>, n: number, k: number): [number[], number[]] {
    log(`Starting findPivots | B: ${B}, S: ${S}, Initial distances: ${d}`);
    const W: Set<number> = new Set<number>(S);
    let W0: number[] = S.slice();
    let Wi: number[] = [];
    for (let i = 1; i <= k; i++) {
        Wi = [];
        const added = new Set<number>();
        for (let u of W0) {
            for (let edge of adj.get(u) || []) {
                const v: number = edge.to;
                const new_d = d[u]! + edge.w;
                const new_depth = depth[u]! + 1;
                const update = new_d < d[v]! || (new_d === d[v] && (new_depth < depth[v]! || (new_depth === depth[v]! && u < Pred[v]!)));
                if (update) {
                    log(`Updating vertex ${v} | New distance: ${new_d}, Previous: ${d[v]}`);
                    d[v] = new_d;
                    depth[v] = new_depth;
                    Pred[v] = u;
                }
                if (d[u]! + edge.w < B && !added.has(v)) {
                    Wi.push(v);
                    added.add(v);
                }
            }
        }
        Wi.forEach(v => W.add(v));
        log(`Iteration ${i} | Wi: ${Wi}, W: ${Array.from(W)}, Distances: ${d}`);
        if (W.size > k * S.length) {
            return [S, Array.from(W)];
        }
        W0 = Wi;
    }
    log(`Ending findPivots | W: ${Array.from(W)}, Final distances: ${d}`);
    return [Array.from(W), Array.from(W)];
}

function baseCase(B: number, S: number[], d: number[], depth: number[], Pred: number[], adj: Map<number, { to: number, w: number }[]>, k: number): [number, number[]] {
    log('Starting baseCase with B:', B, 'S:', S, 'Initial distances:', d);
    const x: number = S[0]!; // Singleton
    const U0: number[] = [x];
    const H: PriorityQueue = new PriorityQueue();
    H.insert(x, d[x]!);
    while (!H.isEmpty() && U0.length < k + 1) {
        const u = H.extractMin()!;
        log('Processing vertex:', u, 'Distance:', d[u]);
        U0.push(u);
        for (let edge of adj.get(u) || []) {
            const v = edge.to;
            const new_d = d[u]! + edge.w;
            if (new_d >= B) continue;
            const new_depth = depth[u]! + 1;
            const update = new_d < d[v]! || (new_d === d[v]! && (new_depth < depth[v]! || (new_depth === depth[v]! && u < Pred[v]!)));
            if (update) {
                log('Updating vertex:', v, 'New distance:', new_d, 'Previous:', d[v]);
                d[v] = new_d;
                depth[v] = new_depth;
                Pred[v] = u;
            }
            if (!H.has(v)) {
                H.insert(v, d[v]!);
            }
        }
    }
    log('Ending baseCase with U0:', U0, 'Final distances:', d);
    if (U0.length <= k) {
        return [B, U0];
    } else {
        const lastD = d[U0[U0.length - 1]!];
        const U = U0.filter(v => d[v]! < lastD!);
        return [lastD!, U];
    }
}

// Step 3: Recursive BMSSP function
function bmssp(l: number, B: number, S: number[], d: number[], depth: number[], Pred: number[], adj: Map<number, { to: number, w: number }[]>, k: number, t: number): [number, number[]] {
    if (l === 0) {
        return baseCase(B, S, d, depth, Pred, adj, k);
    }
    const [P, W] = findPivots(B, S, d, depth, Pred, adj, adj.size, k); // n ~ adj.size, but use graph.n
    const M = Math.pow(2, (l - 1) * t);
    const D = new PartialSortDS(M, B);
    let minP = Infinity;
    for (let x of P) {
        D.insert(x, d[x]!);
        minP = Math.min(minP, d[x]!);
    }
    let B0 = P.length > 0 ? minP : B;
    let U = new Set<number>();
    let i = 0;
    let lastB_prime = B;
    while (U.size < k * Math.pow(2, l * t) && !D.isEmpty()) {
        i++;
        const { x: Bi, S: Si } = D.pull();
        const [B_i_prime, Ui] = bmssp(l - 1, Bi, Si, d, depth, Pred, adj, k, t);
        Ui.forEach(u => U.add(u));
        lastB_prime = B_i_prime;
        const K: { key: number, value: number }[] = [];
        for (let u of Ui) {
            for (let edge of adj.get(u) || []) {
                const v = edge.to;
                const new_d = d[u]! + edge.w;
                const new_depth = depth[u]! + 1;
                const update = new_d < d[v]! || (new_d === d[v]! && (new_depth < depth[v]! || (new_depth === depth[v]! && u < Pred[v]!)));
                if (update) {
                    d[v] = new_d;
                    depth[v] = new_depth;
                    Pred[v] = u;
                }
                if (Bi <= new_d && new_d < B) {
                    D.insert(v, new_d);
                } else if (B_i_prime <= new_d && new_d < Bi) {
                    K.push({ key: v, value: new_d });
                }
            }
        }
        const toPrepend = K;
        for (let x of Si) {
            if (B_i_prime <= d[x]! && d[x]! < Bi) {
                toPrepend.push({ key: x, value: d[x]! });
            }
        }
        D.batchPrepend(toPrepend);
    }
    const B_prime = Math.min(lastB_prime, B);
    const U_arr = Array.from(U);
    W.forEach(x => {
        if (d[x]! < B_prime) U_arr.push(x);
    });
    return [B_prime, U_arr];
}

// Step 4: Main SSSP function
function sssp(graph: Graph, s: number): number[] {
    const n = graph.n;
    const logn = Math.log2(n);
    const k = Math.floor(Math.pow(logn, 1 / 3));
    const t = Math.floor(Math.pow(logn, 2 / 3));
    const l = Math.ceil(logn / t);
    const d = new Array(n).fill(Number.POSITIVE_INFINITY);
    const depth = new Array(n).fill(0);
    const Pred = new Array(n).fill(-1);
    d[s] = 0;
    depth[s] = 1;
    bmssp(l, Number.POSITIVE_INFINITY, [s], d, depth, Pred, graph.adj, k, t);
    return d;
}

// Example of use:
// const graph: Graph = { n: 5, adj: new Map([[0, [{to:1, w:1}, {to:2, w:4}]] , [1, [{to:3, w:2}]], [2, [{to:3, w:1}]] ]) };
// const distances = sssp(graph, 0);
// console.log(distances);

const graph: Graph = {
  n: 4,
  adj: new Map([
    [0, [{to: 1, w: 1}, {to: 2, w: 4}]],
    [1, [{to: 2, w: 2}, {to: 3, w: 5}]],
    [2, [{to: 3, w: 1}]],
  ])
};
const distances = sssp(graph, 0);
console.log(distances);  // Expected: [0, 1, 3, 4] or similar
