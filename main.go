package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
)

// a Graph is with key (nodeIdx) and Node
type Graph map[uint32]*Node

// Node type
type Node struct {
	nodeValue float64
	adjacent  map[uint32]*Node
	weight    map[uint32]float64
}

// Add Node
func (g Graph) AddNode(k uint32) {
	if ptrNode, exist := g[k]; exist {
		err := fmt.Errorf("[Error] Node %v is already exist @%v.", k, ptrNode)
		fmt.Println(err.Error())
	} else {
		g[k] = &Node{nodeValue: 0, adjacent: make(map[uint32]*Node), weight: make(map[uint32]float64)}
	}
}

// Add Edge
func (g Graph) AddEdge(srcIdx, dstIdx uint32, weight float64) {
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
					g.AddNode(i)
				}
			} else {
				for i := uint32(1); i < nNodes+1; i++ {
					g.AddNode(i)
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
				g.AddNode(i)
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

func main() {
	fmt.Println("[Info] Program start.")

	g := Graph{}

	// initialize the graph with SNAP datasets
	initFileName := os.Args[1]
	isWeighted := os.Args[2] == "true"
	//fromZeroIdx := os.Args[3] == "true"
	//g.InitSNAP(initFileName, isWeighted, fromZeroIdx)
	g.InitMatMarket(initFileName, isWeighted)
	//g.Print()
	fmt.Println("[Info] Graph construction done.")

	//// empty list
	//var visitedOrder []uint32
	//
	//// define the Callback function in traversal, here append the search order
	//visitCb := func(i uint32) {
	//	visitedOrder = append(visitedOrder, i)
	//}
	//
	//// start from node 1
	//startInx := uint32(1)
	//DFS(g, startInx, visitCb)
	//fmt.Println(visitedOrder)
	//
	//visitedOrder = []uint32{}
	//
	//BFS(g, startInx, visitCb)
	//fmt.Println(visitedOrder)

	//PageRank(g, 0.85, 0.00000000001)
	//fmt.Println("[Info] PageRank accomplished.")
	maxDist := 0.0
	for v := uint32(1); v < 2; v++ {
		dist, _ := SSSP(g, v)
		for _, val := range dist {
			if val > maxDist {
				maxDist = val
			}
		}
	}
	fmt.Println("[Info] SSSP accomplished.")
	fmt.Println("[Info] MaxDist is:", maxDist)
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
func BFS(g Graph, startIdx uint32, visitCb func(uint32)) {

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

// SSSP
func SSSP(g Graph, startIdx uint32) (map[uint32]float64, map[uint32]bool) {
	visited := map[uint32]bool{}
	dist := map[uint32]float64{}

	visited[startIdx] = true
	dist[startIdx] = 0.0

	for toVisitQueue := []uint32{startIdx}; len(toVisitQueue) > 0; {
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)

		// set the visited mark and deque it from toVisitQueue
		visited[currentNodeIdx] = true
		toVisitQueue = toVisitQueue[1:]

		for v := range currentNode.adjacent {
			ndist := dist[currentNodeIdx] + currentNode.weight[v]
			if !visited[v] {
				toVisitQueue = append(toVisitQueue, v) // bfs
				dist[v] = ndist
			} else if dist[v] > ndist {
				dist[v] = ndist
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
