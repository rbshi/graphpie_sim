package main

import "container/heap"

// SsspPEval
//Note: for the edge-cut scheme, only the worker with startNode should do PEval
func SsspPEval(g Graph, startIdx uint32, nodeRange NodeIDxRange, updateMessage ValueMap, dist ValueMap, visited BoolMap, comPIESend *int, comPIERec *int, workPIE *int) {

	if NodeInRange(startIdx, nodeRange) {
		dist[startIdx] = 0.0
		toVisitQueue := make(PriorityQueue, 0)
		heap.Init(&toVisitQueue)
		heap.Push(&toVisitQueue, &ValueIdxTuple{value: 0.0, nodeIdx: startIdx})
		SsspKernel(g, nodeRange, updateMessage, dist, visited, toVisitQueue, comPIESend, comPIERec, workPIE)
	}
}

func SsspIncEval(g Graph, nodeRange NodeIDxRange, updateMessageMerged ValueMap, updateMessage ValueMap, dist ValueMap, visited BoolMap, comPIESend *int, comPIERec *int, workPIE *int) {

	toVisitQueue := make(PriorityQueue, 0)
	heap.Init(&toVisitQueue)

	//TODO: receive the merged update and write the updated value to dist (DO NOT forget to compare with the local)
	for updateNodeIdx, updateValue := range updateMessageMerged {
		if NodeInRange(updateNodeIdx, nodeRange) {
			distValue, distExist := dist[updateNodeIdx]

			(*comPIERec)++

			if (distExist && updateValue < distValue) || !distExist {
				dist[updateNodeIdx] = updateValue
				heap.Push(&toVisitQueue, &ValueIdxTuple{value: updateValue, nodeIdx: updateNodeIdx})
			}
		}
	}

	SsspKernel(g, nodeRange, updateMessage, dist, visited, toVisitQueue, comPIESend, comPIERec, workPIE)

}

// the overlapped behavior of PEval ^ IncEval: for hardware reuse
func SsspKernel(g Graph, nodeRange NodeIDxRange, updateMessage ValueMap, dist ValueMap, visited BoolMap, toVisitQueue PriorityQueue, comPIESend *int, comPIERec *int, workPIE *int) {
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

			(*workPIE)++

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

					(*comPIESend)++

				}
			}
		}
	}
}

// Sssp
func Sssp(g Graph, startIdx uint32, comSGA *int, workSGA *int) (ValueMap, BoolMap) {
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

			*comSGA++
			*workSGA++

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