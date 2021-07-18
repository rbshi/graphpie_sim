package main

import (
	"fmt"
	"math"
	"os"
)




func main() {
	fmt.Println("[Info] Program start.")

	g := Graph{}

	// initialize the Graph with SNAP datasets
	initFileName := os.Args[1]
	isWeighted := os.Args[2] == "true"
	//fromZeroIdx := os.Args[3] == "true"
	//g.InitSNAP(initFileName, isWeighted, fromZeroIdx)
	g.InitMatMarket(initFileName, isWeighted)
	//g.InitNodeEdgeFile(initFileName, isWeighted)
	//g.Print()
	fmt.Println("[Info] Graph construction done.")

	//resFileName := os.Args[3]
	//res := LoadResult(resFileName)
	//fmt.Println("[Info] Read in the results.")



	// standard

	comSGA := 0
	workSGA := 0

	//SSSP with PIE
	{

		//standard SSSP
 		startNodeIdx := uint32(10000)
		dist, _ := Sssp(g, startNodeIdx, &comSGA, &workSGA)

		fmt.Println("comSGA:", comSGA)
		fmt.Println("workSGA:", workSGA)

		nFrag := 64
		sg := SegmentedPartitioner(g, nFrag)
		fmt.Println("[Info] Number of workers:", len(sg))

		comPIESend := make([]int, nFrag)
		comPIERec := make([]int, nFrag)
		workPIE := make([]int, nFrag)

		//for v, vv := range res {
		//	if dist[v] != vv {
		//		fmt.Println("Different @", v, res[v], vv)
		//	}
		//}

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

			comPIESend[i] = 0
			comPIERec[i] = 0
			workPIE[i] = 0

			SsspPEval(sg[i], startNodeIdx, nodeRange, updateMessage[i], distPtrList[i], visitedPtrList[i], &(comPIESend[i]), &(comPIERec[i]), &(workPIE[i]))
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
				SsspIncEval(sg[i], nodeRange, updateMessageMerged, updateMessage[i], distPtrList[i], visitedPtrList[i], &(comPIESend[i]), &(comPIERec[i]), &(workPIE[i]))
			}

			for i := uint32(0); i < uint32(nFrag); i++ {
				flagNextRound = flagNextRound || len(updateMessage[i]) > 0
			}

			fmt.Println("[Info] Round:", roundIdx)
			roundIdx++
		}

		fmt.Println("[Info] comPIESend:", SliceAccum(comPIESend), "comPIERec:", SliceAccum(comPIERec), "workPIE:", SliceAccum(workPIE))


		for i := uint32(0); i < uint32(nFrag); i++ {
			for v, vv := range distPtrList[i] {
				if dist[v] != vv {
					fmt.Println("Different @", v, dist[v], vv)
				}
			}
		}

		for i := uint32(0); i < uint32(nFrag); i++ {
			fmt.Print(comPIESend[i], "\n")
		}
	}


	// PageRank
	//{
	//	nRound := 10
	//
	//	PageRank(g, 0.85, 0, nRound, &comSGA, &workSGA)
	//
	//	fmt.Println("comSGA:", comSGA)
	//	fmt.Println("workSGA:", workSGA)
	//
	//
	//	nFrag := 64
	//	sg := SegmentedPartitioner(g, nFrag)
	//	// this is the global map for updating message
	//	updateMessage := make([]ValueMap, nFrag)
	//	// PtrList only holds and bypasses the ptr (these two lists are maintained individually by workers)
	//	prPtrList := make([]ValueMap, nFrag)
	//	vnumFrag := uint32(math.Ceil(float64(len(g)) / float64(nFrag)))
	//	damping := 0.85
	//
	//	comPIESend := make([]int, nFrag)
	//	comPIERec := make([]int, nFrag)
	//	workPIE := make([]int, nFrag)
	//
	//	for i := uint32(0); i < uint32(nFrag); i++ {
	//		// FIXME: initialize
	//		prPtrList[i] = ValueMap{}
	//		updateMessage[i] = ValueMap{}
	//		nodeRange := NodeIDxRange{i*vnumFrag + 1, (i + 1) * vnumFrag}
	//		PageRankPEval(sg[i], damping, len(g), nodeRange, updateMessage[i], prPtrList[i], &comPIESend[i], &comPIERec[i], &workPIE[i])
	//	}
	//
	//	//roundIdx := 1
	//	for ii := 0; ii < nRound; ii++ {
	//
	//		//Note: Coordinate updateMessage reduction
	//		updateMessageMerged := ValueMap{}
	//		for v := range g {
	//			for i := uint32(0); i < uint32(nFrag); i++ {
	//				updateValue, updateExist := updateMessage[i][v]
	//				if updateExist {
	//					MapAccum(updateMessageMerged, v, updateValue)
	//				}
	//			}
	//		}
	//		//NOTE: clear the merged message
	//		for i := uint32(0); i < uint32(nFrag); i++ {
	//
	//			comPIESend[i] += len(updateMessage[i])
	//
	//			updateMessage[i] = ValueMap{}
	//		}
	//
	//		//TODO: Assemble (HW)
	//		for i := uint32(0); i < uint32(nFrag); i++ {
	//			nodeRange := NodeIDxRange{i*vnumFrag + 1, (i + 1) * vnumFrag}
	//			for updateNodeIdx, updateValue := range updateMessageMerged {
	//				if NodeInRange(updateNodeIdx, nodeRange) {
	//
	//					comPIERec[i]++
	//
	//					prPtrList[i][updateNodeIdx] += updateValue
	//				}
	//			}
	//		}
	//
	//		fmt.Println("[Info] comPIESend:", SliceAccum(comPIESend), "comPIERec:", SliceAccum(comPIERec), "workPIE:", SliceAccum(workPIE))
	//
	//		//// error calculation
	//		//errSumAllFrag := 0.0
	//		//for v := range g {
	//		//	prValue := prPtrList[int(math.Floor(float64(v-1)/float64(vnumFrag)))][v]
	//		//	errSumAllFrag += math.Abs(g[v].nodeValue - prValue)
	//		//	g[v].nodeValue = prValue
	//		//}
	//		//fmt.Println("[Info] Round:", roundIdx, "AvgError:", errSumAllFrag/float64(len(g)))
	//		//roundIdx++
	//		//
	//		//if errSumAllFrag/float64(len(g)) < 0.000000001 {
	//		//	break
	//		//}
	//
	//		for i := uint32(0); i < uint32(nFrag); i++ {
	//			nodeRange := NodeIDxRange{i*vnumFrag + 1, (i + 1) * vnumFrag}
	//			PageRankKernel(sg[i], damping, len(g), nodeRange, updateMessage[i], prPtrList[i], &comPIESend[i], &comPIERec[i], &workPIE[i])
	//		}
	//	}
	//	fmt.Println("[Info] Accomplished.")
	//	//for i := uint32(0); i < uint32(nFrag); i++ {
	//	//	fmt.Print(workPIE[i], "\n")
	//	//}
	//}

	fmt.Println("[Info] Accomplished.")
}
