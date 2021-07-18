package main

import "math"

// PageRank
func PageRank(g Graph, damping float64, eps float64, nRound int, comSGA *int, workSGA *int) {

	initPR := 1 / float64(len(g))
	for _, v := range g {
		v.nodeValue = initPR
	}

	fixedDump := (1.0 - damping) / float64(len(g))
	//for errAvgSum := 1.0; errAvgSum >= eps; {
	for ii := 0; ii < nRound; ii++ {

		errSum := 0.0

		// calculate the outgoing value of each node
		sendPR := make(ValueMap)
		for i, v := range g {
			if len(v.adjacent) > 0 {
				sendPR[i] = v.nodeValue / float64(len(v.adjacent)) * damping

				*workSGA++

			} else {
				sendPR[i] = 0.0
			}

		}

		for _, v := range g {
			update := 0.0
			for j := range v.adjacent {
				update += sendPR[j]

				*workSGA++
				*comSGA++

			}
			newNodeValue := fixedDump + update
			errSum += math.Abs(newNodeValue - v.nodeValue)
			v.nodeValue = newNodeValue
		}

		// calculate the average error
		//errAvgSum = errSum / float64(len(g))
		//fmt.Println("[Info] PageRank epsilon is converged to: ", errAvgSum)
	}
}

func PageRankPEval(g Graph, damping float64, totalnNode int, nodeRange NodeIDxRange, updateMessage ValueMap, pr ValueMap, comPIESend *int, comPIERec *int, workPIE *int) {

	// average initialize PR
	initPR := 1 / float64(totalnNode)
	for v := range g {
		pr[v] = initPR
		g[v].nodeValue = initPR //FIXME: use nodeValue to keep the old pr for err calculation, NOT used for PR
	}

	PageRankKernel(g, damping, totalnNode, nodeRange, updateMessage, pr, comPIESend, comPIERec, workPIE)

}

func PageRankKernel(g Graph, damping float64, totalnNode int, nodeRange NodeIDxRange, updateMessage ValueMap, pr ValueMap, comPIESend *int, comPIERec *int, workPIE *int) {

	fixedDump := (1.0 - damping) / float64(totalnNode)

	// newPR is for local node
	newPR := ValueMap{}
	for i, v := range g {
		if len(v.adjacent) > 0 {
			// calculate the sendPR
			sendPR := pr[i] / float64(len(v.adjacent)) * damping
			*workPIE++
			// accumulate to newpr tables
			for j := range v.adjacent {
				*workPIE++
				// if local
				if NodeInRange(j, nodeRange) {
					MapAccum(newPR, j, sendPR)
				} else {
					MapAccum(updateMessage, j, sendPR)
				}
			}
		}
	}
	for v := range g {
		pr[v] = newPR[v] + fixedDump
	}
}

