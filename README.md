# cesium heatmap ts

## Project description
```
该库是对开源的kriging 进行的ts版本的编写，没对做其他的改动
原始cesium heatMap的地址 https://github.com/oeo4b/kriging.js
```

## using
```
/**
 * @description: 等值面创建
 * @param {KrigingClass} kriging 克里金对象
 * @param {Array} lngs 经度集合
 * @param {Array} lats 纬度集合
 * @param {Array} siteValue 数值集合
 * @return {*}
 */
const isosurfaces = (kriging: KrigingClass,lngs: Array<number>,lats: Array<number>,siteValue: Array<number>) => {
  const colors = [
    { min: 0, max: 0.1, color: '#FFFFFF' },
    { min: 0.2, max: 10, color: '#A7F290' },
    { min: 11, max: 25, color: '#3CBB3C' },
    { min: 26, max: 50, color: '#61B8FF' },
    { min: 51, max: 100, color: '#0000E1' },
    { min: 101, max: 150, color: '#FA01FA' },
    { min: 151, max: 250, color: '#800040' },
    { min: 251, max: 999, color: '#3F001C' }
  ]

  // 如果存在雨量图则删除雨量图
  const KrigingRain = viewer.entities.getById('KrigingRain')
  viewer.entities.remove(KrigingRain)

  // 绘制面的所有点
  const coords = [] as number[]
  jsonData.geometries[0].coordinates[0].forEach((item:Array<number>) => {
    coords.push(item[0])
    coords.push(item[1])
  })

  // 绘制面的jeojson
  const ex = jsonData.geometries[0].coordinates

  const extent = Cesium.PolygonGeometry.computeRectangle({
    polygonHierarchy: new Cesium.PolygonHierarchy(
      Cesium.Cartesian3.fromDegreesArray(coords)
    )
  }) // 范围（弧度）

  const minx = Cesium.Math.toDegrees(extent.west) // 转换为经纬度
  const miny = Cesium.Math.toDegrees(extent.south)
  const maxx = Cesium.Math.toDegrees(extent.east)
  const maxy = Cesium.Math.toDegrees(extent.north)

  let canvas = null // 画布

  // 1.用克里金训练一个variogram对象
  const variogram = kriging.train(siteValue, lngs, lats, 'spherical', 0, 200)
  // 2.使用刚才的variogram对象使polygons描述的地理位置内的格网元素具备不一样的预测值；
  const grid = kriging.grid(ex, variogram, (maxy - miny) / 1000)

  canvas = document.createElement('canvas')
  canvas.width = 2000
  canvas.height = 2000
  canvas.style.display = 'block'
  // 3.将得到的格网预测值渲染到canvas画布上

  if (grid) {
    kriging.plot(canvas, grid, [minx, maxx], [miny, maxy], colors)
  }

  if (canvas != null) {
    viewer.entities.add({
      id: 'KrigingRain',
      polygon: {
        show: true,
        clampToGround: true,
        hierarchy: {
          positions: Cesium.Cartesian3.fromDegreesArray(coords)
        },
        material: new Cesium.ImageMaterialProperty({
          // 使用贴图的方式将结果贴到面上
          image: canvas,
          transparent: true,
          color: Cesium.Color.WHITE.withAlpha(0.7)
        })
      }
    })
  }
  return viewer.entities.getById('KrigingRain')
}

其中jsonData为GeoJson，根据个人数据进行上述方法的读取插值范围的的代码修改

```

### expectation
```
后续有空再来优化该代码，进行更灵活的使用
```