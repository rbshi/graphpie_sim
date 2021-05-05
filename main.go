package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

// Graph type
//type Graph struct {
//	nodes []*Node
//}

// a Graph is with key (nodeIdx) and Node
type Graph map[int]*Node

// Node type
type Node struct {
	//TODO: should change key to the value of node; add a list to record the edge weights
	key int
	adjacent map[int]*Node
}

// Add Node
func (g Graph) AddNode(k int) {
	if ptrNode, exist := g[k]; exist {
		err := fmt.Errorf("[Error] Node %v is already exist @%v.", k, ptrNode)
		fmt.Println(err.Error())
	} else {
		g[k] = &Node{key: k, adjacent: make(map[int]*Node)}
	}
}

// Add Edge
func (g Graph) AddEdge(srcIdx, dstIdx int) {
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
}

func (g Graph) GetNode(k int) *Node {
	if ptrNode, exist := g[k]; exist {
		return ptrNode
	}
	return nil
}

func (g Graph) InitSNAP(fileName string) {
	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	var nNodes, nEdges int
	for scanner.Scan() {
		//strNodes := strings.Fields(scanner.Text())
		var srcNodeIdx, dstNodeIdx int
		if scanner.Text()[0] != '#' {
			fmt.Sscanf(scanner.Text(), "%v\t%v", &srcNodeIdx, &dstNodeIdx)
			g.AddEdge(srcNodeIdx, dstNodeIdx)
		} else if strings.Contains(scanner.Text(), "Nodes:") {
			fmt.Sscanf(scanner.Text(), "# Nodes: %v Edges: %v", &nNodes, &nEdges)
			// initialize the nodes (bypass the check in AddNode)
			for i := 0; i< nNodes+1; i++ {
				g.AddNode(i)
			}
		}
	}
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
}



func main() {
	fmt.Println("[Info] Program start")

	g := Graph{}

	// initialize the graph with SNAP datasets
	initFileName := os.Args[1]
	g.InitSNAP(initFileName)

	//g.Print()
	fmt.Println("[Info] Graph construction done.")

	// empty list
	var visitedOrder []int

	// define the Callback function in traversal, here append the search order
	visitCb := func(i int) {
		visitedOrder = append(visitedOrder, i)
	}

	// start from node 1
	startInx := 1
	DFS(g, startInx, visitCb)
	fmt.Println(visitedOrder)

	visitedOrder = []int{}

	BFS(g, startInx, visitCb)
	fmt.Println(visitedOrder)
}

func (g Graph) Print() {
	for v, _ := range g {
		fmt.Printf("\nNode %v : ", v)
		for u, _ := range g.GetNode(v).adjacent {
			fmt.Printf(" %v ", u)
		}
	}
}

// DFS
func DFS (g Graph, startIdx int, visitCb func(int)) {
	visited := map[int]bool {}
	visited[startIdx] = true
	for toVisitQueue := []int{startIdx}; len(toVisitQueue) >0; {
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		// task
		visitCb(currentNodeIdx)
		// set the visited mark and deque it from toVisitQueue
		visited[currentNodeIdx] = true
		toVisitQueue = toVisitQueue[1:]

		// preappend the adjacent list of currentNode at the head of the queue (dfs)
		for v, _ := range currentNode.adjacent {
			if !visited[v] {
				toVisitQueue = append([]int{v}, toVisitQueue...)
			}
		}
	}
}


// BFS
func BFS (g Graph, startIdx int, visitCb func(int)) {

	visited := map[int]bool{}

	for toVisitQueue := []int{startIdx}; len(toVisitQueue) >0; {
		// get the currentNodeIdx in Queue
		currentNodeIdx := toVisitQueue[0]
		currentNode := g.GetNode(currentNodeIdx)
		// task
		visitCb(currentNodeIdx)
		// set the visited mark and deque it from toVisitQueue
		visited[currentNodeIdx] = true
		toVisitQueue = toVisitQueue[1:]
		// append the adjacent list of currentNode at the end of queue (bfs)
		for v, _ := range currentNode.adjacent {
			if !visited[v] {
				toVisitQueue = append(toVisitQueue, v)
			}
		}
	}
}
















