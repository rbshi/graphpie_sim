package main



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

