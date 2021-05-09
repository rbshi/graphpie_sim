package main

import (
	"bufio"
	"container/heap"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
)

// Priority queue with heap
type ValueIdxTuple struct {
	value   float64
	nodeIdx uint32
	index   int
}

type PriorityQueue []*ValueIdxTuple

func (pq PriorityQueue) Len() int { return len(pq) }

func (pq PriorityQueue) Less(i, j int) bool {
	// We want Pop to give us the highest, not lowest, priority so we use greater than here.
	return pq[i].value < pq[j].value
}
func (pq PriorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}
func (pq *PriorityQueue) Push(x interface{}) {
	n := pq.Len()
	item := x.(*ValueIdxTuple)
	item.index = n
	*pq = append(*pq, item)
}
func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	old[n-1] = nil  // avoid memory leak
	item.index = -1 // for safety
	*pq = old[0 : n-1]
	return item
}
func (pq *PriorityQueue) update(item *ValueIdxTuple, value float64, nodeIdx uint32) {
	item.value = value
	item.nodeIdx = nodeIdx
	heap.Fix(pq, item.index)
}

// a Graph is with key (nodeIdx) and node
type Graph map[uint32]*Node

// SubGraph is the partitioned from Graph
type SubGraph []Graph

// Node type
type Node struct {
	nodeValue float64
	adjacent  map[uint32]*Node
	weight    ValueMap
}

// visited map
type BoolMap map[uint32]bool
type ValueMap map[uint32]float64

// NodeIdx Range
type NodeIDxRange struct {
	lIdx uint32
	rIdx uint32
}

// Graph partition: SegmentedPartitioner
func SegmentedPartitioner(g Graph, nFrag int) (sg SubGraph) {
	// initialization
	sg = make(SubGraph, nFrag)
	for iFrag := 0; iFrag < nFrag; iFrag++ {
		sg[iFrag] = make(map[uint32]*Node)
	}
	vnumFrag := uint32(math.Ceil(float64(len(g)) / float64(nFrag)))
	for srcNodeIdx, srcNode := range g {
		idxFrag := uint32(math.Floor(float64(srcNodeIdx-1) / float64(vnumFrag)))
		sg[idxFrag][srcNodeIdx] = srcNode
	}
	return sg
}

// Add Node
func (g Graph) AddNode(k uint32, v float64) {
	if ptrNode, exist := g[k]; exist {
		err := fmt.Errorf("[Error] Node %v is already exist @%v.", k, ptrNode)
		fmt.Println(err.Error())
	} else {
		g[k] = &Node{nodeValue: v, adjacent: make(map[uint32]*Node), weight: make(ValueMap)}
	}
}

// Add Edge
func (g Graph) AddEdge(srcIdx, dstIdx uint32, weight float64) {
	// get the Node
	srcNode := g.GetNode(srcIdx)
	dstNode := g.GetNode(dstIdx)
	// check if the node exist
	if srcNode == nil || dstNode == nil {
		err := fmt.Errorf("[Error] Edge invalid: (%v -> %v).", srcIdx, dstIdx)
		fmt.Println(err.Error())
	} else if _, exist := srcNode.adjacent[dstIdx]; exist {
		err := fmt.Errorf("[Error] Edge already exists: (%v -> %v).", srcIdx, dstIdx)
		fmt.Println(err.Error())
	}
	// add the edge
	srcNode.adjacent[dstIdx] = dstNode
	srcNode.weight[dstIdx] = weight
}

func (g Graph) GetNode(k uint32) *Node {
	if ptrNode, exist := g[k]; exist {
		return ptrNode
	}
	return nil
}

func (g Graph) InitSNAP(fileName string, isWeighted bool, fromZeroIdx bool) {
	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	var nNodes, nEdges uint32
	for scanner.Scan() {
		//strNodes := strings.Fields(scanner.Text())
		var srcNodeIdx, dstNodeIdx uint32
		edgeWeight := 1.0
		if scanner.Text()[0] != '#' {
			if isWeighted {
				fmt.Sscanf(scanner.Text(), "%v\t%v\t%v", &srcNodeIdx, &dstNodeIdx, &edgeWeight)
			} else {
				fmt.Sscanf(scanner.Text(), "%v\t%v", &srcNodeIdx, &dstNodeIdx)
			}
			g.AddEdge(srcNodeIdx, dstNodeIdx, edgeWeight)
		} else if strings.Contains(scanner.Text(), "Nodes:") {
			fmt.Sscanf(scanner.Text(), "# Nodes: %v Edges: %v", &nNodes, &nEdges)
			// initialize the nodes (bypass the check in AddNode)
			//TODO: the algorithm codes are with NodeIdx from 1 (with Florida format)
			if fromZeroIdx {
				for i := uint32(0); i < nNodes; i++ {
					g.AddNode(i, 0)
				}
			} else {
				for i := uint32(1); i < nNodes+1; i++ {
					g.AddNode(i, 0)
				}
			}
		}
	}
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
}

// For Florida sparse matrix collection with Matrix Market format
func (g Graph) InitMatMarket(fileName string, isWeighted bool) {
	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	var nNodes, nEdges uint32
	for scanner.Scan() {
		// bypass the comments
		if scanner.Text()[0] == '%' {
			continue
		} else {
			fmt.Sscanf(scanner.Text(), "%v %v %v", &nNodes, &nNodes, &nEdges)
			// add nodes
			for i := uint32(1); i <= nNodes; i++ {
				g.AddNode(i, 0)
			}
			break
		}
	}

	var srcNodeIdx, dstNodeIdx uint32
	edgeWeight := 1.0
	for scanner.Scan() {
		//strNodes := strings.Fields(scanner.Text())
		if isWeighted {
			fmt.Sscanf(scanner.Text(), "%v\t%v\t%v", &srcNodeIdx, &dstNodeIdx, &edgeWeight)
		} else {
			fmt.Sscanf(scanner.Text(), "%v\t%v", &srcNodeIdx, &dstNodeIdx)
		}
		g.AddEdge(srcNodeIdx, dstNodeIdx, edgeWeight)
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
}

func LoadResult(fileName string) (res ValueMap) {
	res = ValueMap{}
	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	nodeValue := 1.0
	nodeIdx := uint32(1)

	for scanner.Scan() {
		// if the node result value is `infinity`, bypass adding it into result map
		if strings.Contains(scanner.Text(), "infinity") {
			continue
		}
		fmt.Sscanf(scanner.Text(), "%v %v", &nodeIdx, &nodeValue)
		res[nodeIdx] = nodeValue
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	return res
}

// For Florida sparse matrix collection with Matrix Market format
func (g Graph) InitNodeEdgeFile(fileName string, isWeighted bool) {
	fileNode, err := os.Open(fileName + ".v")
	if err != nil {
		log.Fatal(err)
	}
	defer fileNode.Close()

	nodeValue := 1.0
	nodeIdx := uint32(1)
	scanner := bufio.NewScanner(fileNode)
	for scanner.Scan() {
		fmt.Sscanf(scanner.Text(), "%v %v", &nodeIdx, &nodeValue)
		g.AddNode(nodeIdx, nodeValue)
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	fileEdge, err := os.Open(fileName + ".e")
	if err != nil {
		log.Fatal(err)
	}
	defer fileEdge.Close()

	var srcNodeIdx, dstNodeIdx uint32
	edgeWeight := 1.0
	scanner = bufio.NewScanner(fileEdge)
	for scanner.Scan() {
		//strNodes := strings.Fields(scanner.Text())
		if isWeighted {
			fmt.Sscanf(scanner.Text(), "%v %v %v", &srcNodeIdx, &dstNodeIdx, &edgeWeight)
		} else {
			fmt.Sscanf(scanner.Text(), "%v\t%v", &srcNodeIdx, &dstNodeIdx)
		}
		g.AddEdge(srcNodeIdx, dstNodeIdx, edgeWeight)
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
}

func main() {
	fmt.Println("[Info] Program start.")

	g := Graph{}

	// initialize the Graph with SNAP datasets
	initFileName := os.Args[1]
	isWeighted := os.Args[2] == "true"
	//fromZeroIdx := os.Args[3] == "true"
	//g.InitSNAP(initFileName, isWeighted, fromZeroIdx)
	//g.InitMatMarket(initFileName, isWeighted)
	g.InitNodeEdgeFile(initFileName, isWeighted)
	//g.Print()
	fmt.Println("[Info] Graph construction done.")

	resFileName := os.Args[3]
	res := LoadResult(resFileName)
	fmt.Println("[Info] Read in the results.")

	nFrag := 16
	sg := SegmentedPartitioner(g, nFrag)
	fmt.Println("[Info] Number of workers:", len(sg))

	startNodeIdx := uint32(6)
	// standard
	dist, _ := Sssp(g, startNodeIdx)

	for v, vv := range res {
		if dist[v] != vv {
			fmt.Println("Different @", v, res[v], vv)
		}
	}

	// this is the global map for updating message
	updateMessage := make([]ValueMap, nFrag)

	// PtrList only holds and bypasses the ptr (these two lists are maintained individually by workers)
	distPtrList := make([]ValueMap, nFrag)
	visitedPtrList := make([]BoolMap, nFrag)

	vnumFrag := uint32(math.Ceil(float64(len(g)) / float64(nFrag)))

	for i := uint32(0); i < uint32(nFrag); i++ {
		// FIXME: initialize
		distPtrList[i] = ValueMap{}
		visitedPtrList[i] = BoolMap{}
		updateMessage[i] = ValueMap{}

		nodeRange := NodeIDxRange{i*vnumFrag + 1, (i + 1) * vnumFrag}

		SsspPEval(sg[i], startNodeIdx, nodeRange, updateMessage[i], distPtrList[i], visitedPtrList[i])
	}

	roundIdx := 1
	flagNextRound := true
	for flagNextRound {

		flagNextRound = false

		//Note: Coordinate updateMessage reduction
		updateMessageMerged := ValueMap{}
		for v := range g {
			m := 100000000.0 //FIXME: what is inf?
			for i := uint32(0); i < uint32(nFrag); i++ {
				updateValue, updateExist := updateMessage[i][v]
				if updateExist && updateValue < m {
					m = updateValue
					updateMessageMerged[v] = m
				}
			}
		}

		//NOTE: clear the merged message
		for i := uint32(0); i < uint32(nFrag); i++ {
			updateMessage[i] = ValueMap{}
		}

		for i := uint32(0); i < uint32(nFrag); i++ {
			// NOTE: should clear the visit in each IncEval? LibGrape did that
			visitedPtrList[i] = BoolMap{}
			nodeRange := NodeIDxRange{i*vnumFrag + 1, (i + 1) * vnumFrag}
			SsspIncEval(sg[i], nodeRange, updateMessageMerged, updateMessage[i], distPtrList[i], visitedPtrList[i])
		}

		for i := uint32(0); i < uint32(nFrag); i++ {
			flagNextRound = flagNextRound || len(updateMessage[i]) > 0
		}

		fmt.Println("[Info] Round:", roundIdx)
		roundIdx++
	}

	for i := uint32(0); i < uint32(nFrag); i++ {
		for v, vv := range distPtrList[i] {
			if dist[v] != vv {
				fmt.Println("Different @", v, dist[v], vv)
			}
		}
	}

	fmt.Println("[Info] Accomplished.")
}

func (g Graph) Print() {
	for v := range g {
		fmt.Printf("\nNode %v : ", v)
		for u := range g.GetNode(v).adjacent {
			fmt.Printf(" %v ", u)
		}
	}
}

// DFS
func DFS(g Graph, startIdx uint32, visitCb func(uint32)) {
	visited := BoolMap{}
	visited[startIdx] = true
	for toVisitQueue := []uint32{startIdx}; len(toVisitQueue) > 0; {
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		// task
		visitCb(currentNodeIdx)
		// set the visited mark and deque it from toVisitQueue
		visited[currentNodeIdx] = true
		toVisitQueue = toVisitQueue[1:]

		// pre-append the adjacent list of currentNode at the head of the queue (dfs)
		for v := range currentNode.adjacent {
			if !visited[v] {
				toVisitQueue = append([]uint32{v}, toVisitQueue...)
			}
		}
	}
}

// BFS
func BFS(g Graph, startIdx uint32, visitCb func(uint32)) {

	visited := BoolMap{}

	for toVisitQueue := []uint32{startIdx}; len(toVisitQueue) > 0; {
		// get the currentNodeIdx in Queue
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		// task
		visitCb(currentNodeIdx)
		// set the visited mark and deque it from toVisitQueue
		visited[currentNodeIdx] = true
		toVisitQueue = toVisitQueue[1:]
		// append the adjacent list of currentNode at the end of queue (bfs)
		for v := range currentNode.adjacent {
			if !visited[v] {
				toVisitQueue = append(toVisitQueue, v)
			}
		}
	}
}

func NodeInRange(nodeIdx uint32, nodeRange NodeIDxRange) bool {
	return nodeIdx >= nodeRange.lIdx && nodeIdx <= nodeRange.rIdx
}

func MinFloat64(v []float64) float64 {
	m := 100000000.0 //FIXME: what is inf?
	for i, e := range v {
		if i == 0 || e < m {
			m = e
		}
	}
	return m
}

// SsspPEval
//Note: for the edge-cut scheme, only the worker with startNode should do PEval
func SsspPEval(g Graph, startIdx uint32, nodeRange NodeIDxRange, updateMessage ValueMap, dist ValueMap, visited BoolMap) {

	if NodeInRange(startIdx, nodeRange) {
		dist[startIdx] = 0.0
		toVisitQueue := make(PriorityQueue, 0)
		heap.Init(&toVisitQueue)
		heap.Push(&toVisitQueue, &ValueIdxTuple{value: 0.0, nodeIdx: startIdx})
		SsspKernel(g, nodeRange, updateMessage, dist, visited, toVisitQueue)
	}
}

func SsspIncEval(g Graph, nodeRange NodeIDxRange, updateMessageMerged ValueMap, updateMessage ValueMap, dist ValueMap, visited BoolMap) {

	toVisitQueue := make(PriorityQueue, 0)
	heap.Init(&toVisitQueue)

	//TODO: receive the merged update and write the updated value to dist (DO NOT forget to compare with the local)
	for updateNodeIdx, updateValue := range updateMessageMerged {
		if NodeInRange(updateNodeIdx, nodeRange) {
			distValue, distExist := dist[updateNodeIdx]
			if (distExist && updateValue < distValue) || !distExist {
				dist[updateNodeIdx] = updateValue
				heap.Push(&toVisitQueue, &ValueIdxTuple{value: updateValue, nodeIdx: updateNodeIdx})
			}
		}
	}

	SsspKernel(g, nodeRange, updateMessage, dist, visited, toVisitQueue)

}

// the overlapped behavior of PEval ^ IncEval: for hardware reuse
func SsspKernel(g Graph, nodeRange NodeIDxRange, updateMessage ValueMap, dist ValueMap, visited BoolMap, toVisitQueue PriorityQueue) {
	for len(toVisitQueue) > 0 {
		currentNodeTuple := heap.Pop(&toVisitQueue).(*ValueIdxTuple)
		currentNode := g.GetNode(currentNodeTuple.nodeIdx)

		// Dijkstra
		if visited[currentNodeTuple.nodeIdx] {
			continue
		}
		visited[currentNodeTuple.nodeIdx] = true

		for v := range currentNode.adjacent {
			ndist := dist[currentNodeTuple.nodeIdx] + currentNode.weight[v]
			// the update will be consumed locally (preferable)
			if NodeInRange(v, nodeRange) {
				_, exist := dist[v]
				if !exist || dist[v] > ndist {
					dist[v] = ndist
					heap.Push(&toVisitQueue, &ValueIdxTuple{value: ndist, nodeIdx: v})
				}
			} else {
				//FIXME: sendout <updated v, ndist>, but compare the update to the same v in the message list (is that possible on HW?)
				updateValue, updateExist := updateMessage[v]
				if !updateExist || updateValue > ndist {
					updateMessage[v] = ndist
				}
			}
		}
	}
}

// Sssp
func Sssp(g Graph, startIdx uint32) (ValueMap, BoolMap) {
	visited := BoolMap{}
	dist := ValueMap{}

	dist[startIdx] = 0.0

	toVisitQueue := make(PriorityQueue, 0)
	heap.Init(&toVisitQueue)
	heap.Push(&toVisitQueue, &ValueIdxTuple{value: 0.0, nodeIdx: startIdx})

	for len(toVisitQueue) > 0 {
		currentNodeTuple := heap.Pop(&toVisitQueue).(*ValueIdxTuple)
		currentNode := g.GetNode(currentNodeTuple.nodeIdx)

		// Dijkstra
		if visited[currentNodeTuple.nodeIdx] {
			continue
		}
		visited[currentNodeTuple.nodeIdx] = true

		for v := range currentNode.adjacent {
			ndist := dist[currentNodeTuple.nodeIdx] + currentNode.weight[v]
			_, exist := dist[v]
			if !exist || dist[v] > ndist {
				dist[v] = ndist
				heap.Push(&toVisitQueue, &ValueIdxTuple{value: ndist, nodeIdx: v})
			}
		}
	}
	return dist, visited
}

// PageRank
func PageRank(g Graph, damping float64, eps float64) {
	sumEdgeWeight := 0.0
	// initialization (PR(node) = 1 / sum(edgeWeight))
	for _, v := range g {
		for _, e := range v.weight {
			sumEdgeWeight += e
		}
	}

	//TODO: make the PR trans in weighted version
	initPR := 1 / float64(len(g))
	for _, v := range g {
		v.nodeValue = initPR
	}

	fixedDump := (1.0 - damping) / float64(len(g))
	for errAvgSum := 1.0; errAvgSum >= eps; {
		errSum := 0.0

		// calculate the outgoing value of each node
		sendNodeValue := make(ValueMap)
		for i, v := range g {
			if len(v.adjacent) > 0 {
				sendNodeValue[i] = v.nodeValue / float64(len(v.adjacent)) * damping
			} else {
				sendNodeValue[i] = 0.0
			}

		}

		for _, v := range g {
			update := 0.0
			for j := range v.adjacent {
				update += sendNodeValue[j]
			}
			newNodeValue := fixedDump + update
			errSum += math.Abs(newNodeValue - v.nodeValue)
			v.nodeValue = newNodeValue
		}

		// calculate the average error
		errAvgSum = errSum / float64(len(g))
		fmt.Println("[Info] PageRank epsilon is converged to: ", errAvgSum)
	}
}
