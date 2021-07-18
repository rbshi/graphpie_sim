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


func NodeInRange(nodeIdx uint32, nodeRange NodeIDxRange) bool {
	return nodeIdx >= nodeRange.lIdx && nodeIdx <= nodeRange.rIdx
}

func (g Graph) Print() {
	for v := range g {
		fmt.Printf("\nNode %v : ", v)
		for u := range g.GetNode(v).adjacent {
			fmt.Printf(" %v ", u)
		}
	}
}

func SliceAccum(s []int) (sum int){
	sum = 0
	for i := range s {
		sum += s[i]
	}
	return sum
}



func MapAccum(valueMap ValueMap, dst uint32, accumValue float64) {
	oriValue, exist := valueMap[dst]
	if exist {
		valueMap[dst] = oriValue + accumValue
	} else {
		valueMap[dst] = accumValue
	}
}