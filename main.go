package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
)

// a graph is with key (nodeIdx) and node
type graph map[uint32]*node

// subGraph is the partitioned from graph
type subGraph []graph

// node type
type node struct {
	nodeValue float64
	adjacent  map[uint32]*node
	weight    map[uint32]float64
}

// visited map
type visitedMap map[uint32]bool
type valueMap map[uint32]float64

// Graph partition: SegmentedPartitioner
func SegmentedPartitioner(g graph, nFrag int) (sg subGraph) {
	// initialization
	sg = make(subGraph, nFrag)
	for iFrag := 0; iFrag < nFrag; iFrag++ {
		sg[iFrag] = make(map[uint32]*node)
	}
	vnumFrag := uint32(math.Ceil(float64(len(g)) / float64(nFrag)))
	for srcNodeIdx, srcNode := range g {
		idxFrag := uint32(math.Floor(float64(srcNodeIdx-1) / float64(vnumFrag)))
		sg[idxFrag][srcNodeIdx] = srcNode
	}
	return sg
}

// Add node
func (g graph) AddNode(k uint32, v float64) {
	if ptrNode, exist := g[k]; exist {
		err := fmt.Errorf("[Error] Node %v is already exist @%v.", k, ptrNode)
		fmt.Println(err.Error())
	} else {
		g[k] = &node{nodeValue: v, adjacent: make(map[uint32]*node), weight: make(map[uint32]float64)}
	}
}

// Add Edge
func (g graph) AddEdge(srcIdx, dstIdx uint32, weight float64) {
	// get the node
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

func (g graph) GetNode(k uint32) *node {
	if ptrNode, exist := g[k]; exist {
		return ptrNode
	}
	return nil
}

func (g graph) InitSNAP(fileName string, isWeighted bool, fromZeroIdx bool) {
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
func (g graph) InitMatMarket(fileName string, isWeighted bool) {
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

func LoadResult(fileName string) (res valueMap) {
	res = valueMap{}
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
func (g graph) InitNodeEdgeFile(fileName string, isWeighted bool) {
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

	g := graph{}

	// initialize the graph with SNAP datasets
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

	nFrag := 2048
	sg := SegmentedPartitioner(g, nFrag)
	fmt.Println("[Info] Number of workers:", len(sg))

	//// empty list
	//var visitedOrder []uint32
	//
	//// define the Callback function in traversal, here append the search order
	//visitCb := func(i uint32) {
	//	visitedOrder = append(visitedOrder, i)
	//}

	startNodeIdx := uint32(6)
	// standard
	dist, _ := Sssp(g, startNodeIdx)

	for v, vv := range res {
		if dist[v] != vv {
			fmt.Println("Different @", v, res[v], vv)
		}
	}

	// this is the global map for updating message
	// FIXME: []float64 is the the message (ndist) in SSSP
	updateMessage := make([]map[uint32]float64, nFrag)

	// PtrList only holds and bypasses the ptr (these two lists are mantained individually by workers)
	distPtrList := make([]map[uint32]float64, nFrag)
	visitedPtrList := make([]map[uint32]bool, nFrag)

	vnumFrag := uint32(math.Ceil(float64(len(g)) / float64(nFrag)))

	for i := uint32(0); i < uint32(nFrag); i++ {
		// FIXME: initialize
		distPtrList[i] = valueMap{}
		visitedPtrList[i] = visitedMap{}
		updateMessage[i] = valueMap{}
		SsspPEval(sg[i], startNodeIdx, i*vnumFrag+1, (i+1)*vnumFrag, updateMessage[i], distPtrList[i], visitedPtrList[i])
	}

	flagNextRound := true
	for flagNextRound {

		flagNextRound = false

		//FIXME: Coordinate updateMessage reduction
		updateMessageMerged := valueMap{}
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

		//FIXME: clear message
		for i := uint32(0); i < uint32(nFrag); i++ {
			updateMessage[i] = valueMap{}
		}

		for i := uint32(0); i < uint32(nFrag); i++ {
			// FIXME: should clear the visit in each IncEval? LibGrape did that
			visitedPtrList[i] = visitedMap{}
			SsspIncEval(sg[i], startNodeIdx, i*vnumFrag+1, (i+1)*vnumFrag, updateMessageMerged, updateMessage[i], distPtrList[i], visitedPtrList[i])
		}

		for i := uint32(0); i < uint32(nFrag); i++ {
			flagNextRound = flagNextRound || len(updateMessage[i]) > 0
		}

		fmt.Println("[Info] One round.")
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

func (g graph) Print() {
	for v := range g {
		fmt.Printf("\nNode %v : ", v)
		for u := range g.GetNode(v).adjacent {
			fmt.Printf(" %v ", u)
		}
	}
}

// DFS
func DFS(g graph, startIdx uint32, visitCb func(uint32)) {
	visited := map[uint32]bool{}
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
func BFS(g graph, startIdx uint32, visitCb func(uint32)) {

	visited := map[uint32]bool{}

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

func NodeInRange(nodeIdx uint32, lRangeVIdx uint32, rRangeVIdx uint32) bool {
	return nodeIdx >= lRangeVIdx && nodeIdx <= rRangeVIdx
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
//TODO: for the vertex partition scheme, only the worker with startNode should do PEval
func SsspPEval(g graph, startIdx uint32, lRangeVIdx uint32, rRangeVIdx uint32, updateMessage valueMap, dist valueMap, visited visitedMap) {

	dist[startIdx] = 0.0

	if NodeInRange(startIdx, lRangeVIdx, rRangeVIdx) {
		for toVisitQueue := []uint32{startIdx}; len(toVisitQueue) > 0; {
			currentNodeIdx := toVisitQueue[0]
			currentNode := g.GetNode(currentNodeIdx)
			toVisitQueue = toVisitQueue[1:]

			//if visited[currentNodeIdx] {
			//	continue
			//}

			visited[currentNodeIdx] = true

			for v := range currentNode.adjacent {
				ndist := dist[currentNodeIdx] + currentNode.weight[v]
				// the update will be consumed locally (preferable)
				if NodeInRange(v, lRangeVIdx, rRangeVIdx) {
					_, exist := dist[v]
					if !exist || dist[v] > ndist {
						dist[v] = ndist
						//FIXME: here still will be multi insert of the same node? but will be only visited ONCE as the src
						toVisitQueue = append(toVisitQueue, v)
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
}

func SsspIncEval(g graph, startIdx uint32, lRangeVIdx uint32, rRangeVIdx uint32, updateMessageMerged valueMap, updateMessage valueMap, dist valueMap, visited visitedMap) {

	toVisitQueue := make([]uint32, 0)

	//TODO: receive the merged update and write the updated value to dist (DO NOT forget to compare with the local)
	for i := lRangeVIdx; i <= rRangeVIdx; i++ {
		updateValue, updateExist := updateMessageMerged[i]
		distValue, distExist := dist[i]
		//FIXME
		if distExist && updateExist && updateValue < distValue {
			dist[i] = updateValue
			toVisitQueue = append(toVisitQueue, i)
		} else if !distExist && updateExist {
			dist[i] = updateValue
			toVisitQueue = append(toVisitQueue, i)
		}
	}

	for len(toVisitQueue) > 0 {
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		toVisitQueue = toVisitQueue[1:]

		//if visited[currentNodeIdx] {
		//	continue
		//}

		visited[currentNodeIdx] = true

		for v := range currentNode.adjacent {
			ndist := dist[currentNodeIdx] + currentNode.weight[v]
			// the update will be consumed locally (preferable)
			if NodeInRange(v, lRangeVIdx, rRangeVIdx) {
				_, exist := dist[v]
				if !exist || dist[v] > ndist {
					dist[v] = ndist
					//FIXME: here still will be multi insert of the same node? but will be only visited ONCE as the src
					toVisitQueue = append(toVisitQueue, v)
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
func Sssp(g graph, startIdx uint32) (valueMap, visitedMap) {
	visited := visitedMap{}
	dist := valueMap{}

	dist[startIdx] = 0.0

	for toVisitQueue := []uint32{startIdx}; len(toVisitQueue) > 0; {
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		toVisitQueue = toVisitQueue[1:]

		//if visited[currentNodeIdx] {
		//	continue
		//}

		visited[currentNodeIdx] = true

		for v := range currentNode.adjacent {
			//FIXME: only for v without being visited!
			ndist := dist[currentNodeIdx] + currentNode.weight[v]
			_, exist := dist[v]
			if !exist || dist[v] > ndist {
				dist[v] = ndist
				//FIXME: here still will be multi insert of the same node? but will be only visited ONCE as the src
				toVisitQueue = append(toVisitQueue, v)
			}
		}
	}
	return dist, visited
}

// PageRank
func PageRank(g graph, damping float64, eps float64) {
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
		sendNodeValue := make(map[uint32]float64)
		for i, v := range g {
			if len(v.adjacent) > 0 {
				sendNodeValue[i] = v.nodeValue / float64(len(v.adjacent)) * damping
			} else {
				sendNodeValue[i] = 0.0
			}

		}

		for _, v := range g {
			update := 0.0
			for j, _ := range v.adjacent {
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
