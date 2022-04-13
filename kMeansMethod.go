package main

import (
	"fmt"
	"math"
)

func calcDistance(point1 []float64, point2 []float64, order int) float64 {
	var squareSum float64
	for i := 0; i < order; i++ {
		squareSum += (point1[i] - point2[i]) * (point1[i] - point2[i])
	}
	return math.Sqrt(squareSum)
}

// 计算源数据一个点到所有重心点的距离
func calcSingleDataDistances(data []float64, order int, gravityPoints [][]float64, setNum int) []float64 {
	arrayDistance := make([]float64, setNum)
	for i := 0; i < setNum; i++ {
		arrayDistance[i] = calcDistance(data, gravityPoints[i], order)
	}
	return arrayDistance
}

func calcAllDataDistances(datas [][]float64, dataNum int, order int, gravityPoints [][]float64, setNum int) [][]float64 {
	allDistances := make([][]float64, dataNum)
	for i := 0; i < dataNum; i++ {
		allDistances[i] = calcSingleDataDistances(datas[i], order, gravityPoints, setNum)
	}
	return allDistances
}

// 计算源数据所有点属于第几个聚类
func calcSetIndex(distances [][]float64, dataNum int) []int {
	indexSets := make([]int, dataNum) // 第几个数据就代表源数据第几个数据的聚类下标
	for i := 0; i < dataNum; i++ {
		var minDis float64
		minDis = distances[i][0]
		indexSets[i] = 0
		for j := 1; j < len(distances[i]); j++ {
			if distances[i][j] < minDis {
				minDis = distances[i][j]
				indexSets[i] = j
			}
		}
	}
	return indexSets
}

func calcNewGravityPoints(datas [][]float64, dataNum int, order int, indexSet []int, setNum int) [][]float64 {
	// 源数据根据现有重心分组
	newdatas := make([][][]float64, setNum)
	for i := 0; i < dataNum; i++ {
		newdatas[indexSet[i]] = append(newdatas[indexSet[i]], datas[i])
	}
	newGravityPoints := make([][]float64, setNum)
	for i := 0; i < setNum; i++ {
		gravityPoint := make([]float64, order)
		for j := 0; j < len(newdatas[i]); j++ {
			for k := 0; k < order; k++ {
				gravityPoint[k] += newdatas[i][j][k]
			}
		}
		for l := 0; l < order; l++ {
			gravityPoint[l] = gravityPoint[l] / float64(len(newdatas[l]))
		}
		newGravityPoints[i] = gravityPoint
	}
	return newGravityPoints
}

// 判断新的重心点是否收敛
func calcGravityError(oldGravityPoints, newGravityPoints [][]float64, order int, setNum int, zero float64) bool {
	errs := make([]float64, setNum)
	for i := 0; i < setNum; i++ {
		err := 0.0
		for j := 0; j < order; j++ {
			err += (oldGravityPoints[i][j] - newGravityPoints[i][j]) * (oldGravityPoints[i][j] - newGravityPoints[i][j])
		}
		errs[i] = err
	}

	for i := 0; i < setNum; i++ {
		if errs[i] > zero {
			return false // 任意一个重心点未收敛返回false
		}
	}
	return true
}

//
func kMeans(datas [][]float64, dataNum int, order int, initGravityPoints [][]float64, setNum int) [][]float64 {
	const zero = 1e3
	var isEnd = false
	for {
		if isEnd {
			break
		}
		// 根据初始重心点计算所有点到各重心的距离 返回值 dataNum * setNum 矩阵
		allDistances := calcAllDataDistances(datas, dataNum, order, initGravityPoints, setNum)
		// 计算所有源数据属于哪一个聚类
		dataIndex := calcSetIndex(allDistances, dataNum)
		// 计算新的重心点
		newGravityPoints := calcNewGravityPoints(datas, dataNum, order, dataIndex, setNum)
		isEnd = calcGravityError(initGravityPoints, newGravityPoints, order, setNum, zero)
		initGravityPoints = newGravityPoints
	}
	return initGravityPoints
}

// an example
func main() {
	datas := [][]float64{
		{1.1, 1.1}, {5.2, 5.3}, {1.8, 0.1}, {4.9, 4.3},
		{1.9, 1.9}, {4.9, 5.2}, {1.2, 1.5}, {5.2, 5.8},
		{1.2, 0.9}, {5.1, 4.8}, {0.9, 1.1}, {6.1, 5.9},
	}
	initGravityPoints := [][]float64{
		{1, 1},
		{5, 5},
	}
	result := kMeans(datas, 12, 2, initGravityPoints, 2)
	fmt.Println(result)
}
