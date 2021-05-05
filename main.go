package main

import (
	"fmt"
)

// Graph type
type Graph struct {
	vertices []*Vertex

}

// Vertex type
type Vertex struct {
	key int
	adjacent []*Vertex
}


// Add Vertex
func (g *Graph) AddVertex(k int) {
	if Contains(g.vertices, k) {
		err := fmt.Errorf("Vertex %v is already exist.", k)
		fmt.Println(err.Error())
	} else {
		g.vertices = append(g.vertices, &Vertex{key: k})
	}
}

// Add Edge
func (g *Graph) AddEdge(src, dst int) {
	// get the vertex
	srcVertex := g.GetVertex(src)
	dstVertex := g.GetVertex(dst)

	// check if the node exist
	if srcVertex == nil || dstVertex == nil {
		err := fmt.Errorf("Edge invalid: (%v -> %v).", src, dst)
		fmt.Println(err.Error())
	} else if Contains(srcVertex.adjacent, dst) {
		err := fmt.Errorf("Edge already exists: (%v -> %v).", src, dst)
		fmt.Println(err.Error())
	}
	// add edge
	srcVertex.adjacent = append(srcVertex.adjacent, dstVertex)
}

func (g *Graph) GetVertex(k int) *Vertex {
	for i, v := range g.vertices {
		if v.key == k {
			return g.vertices[i]
		}
	}
	return nil
}

func Contains(s []*Vertex, k int) bool {
	for _, v := range s{
		if k == v.key {
			return true
		}
	}
	return false
}



func main() {
	fmt.Println("hello world!")

	g := &Graph{}




	for i := 1; i < 11; i++ {
		g.AddVertex(i)
	}

	g.AddEdge(1, 9)
	g.AddEdge(1, 5)
	g.AddEdge(1, 2)
	g.AddEdge(2, 2)
	g.AddEdge(3, 4)
	g.AddEdge(5, 6)
	g.AddEdge(5, 8)
	g.AddEdge(6, 7)
	g.AddEdge(9, 10)

	g.Print()
	fmt.Println("[Done] Graph construction.")


	// empty list
	var visitedOrder []int

	// define the Callback function in traversal, here append the search order
	visitCb := func(i int) {
		visitedOrder = append(visitedOrder, i)
	}

	// start from vertex 1
	startInx := 1
	DFS(g, startInx, visitCb)
	fmt.Println(visitedOrder)

	visitedOrder = []int{}

	BFS(g, startInx, visitCb)
	fmt.Println(visitedOrder)


}

func (g *Graph) Print() {
	for _, v := range g.vertices {
		fmt.Printf("\nVertex %v : ", v.key)
		for _, v := range v.adjacent {
			fmt.Printf(" %v ", v.key)
		}
	}
}


// DFS
func DFS (g *Graph, startIdx int, visitCb func(int)) {
	visited := map[int]bool {}

	startVertex := g.GetVertex(startIdx)
	if startVertex == nil {
		err := fmt.Errorf("Vertex %v does not exist.", startIdx)
		fmt.Println(err.Error())
	} else {

		visited[startVertex.key] = true
		visitCb(startVertex.key)

		for _, v := range startVertex.adjacent {
			if visited[v.key] {
				continue
			}
			DFS(g, v.key, visitCb)
		}
	}
}


// BFS
func BFS (g *Graph, startIdx int, visitCb func(int)) {

	visited := map[int]bool{}

	for toVisitQueue := []int{startIdx}; len(toVisitQueue) >0; {
		// get the currentVertex from the key in queue
		currentVertex := g.GetVertex(toVisitQueue[0])
		visitCb(currentVertex.key)
		// set the visited mark and deque the key from toVisitQueue
		visited[currentVertex.key] = true
		toVisitQueue = toVisitQueue[1:]

		for _, v := range currentVertex.adjacent {
			if !visited[v.key] {
				toVisitQueue = append(toVisitQueue, v.key)
			}
		}
	}
}
















