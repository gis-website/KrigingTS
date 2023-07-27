/* eslint no-extend-native: ["error", { "exceptions": ["Array"] }] */
// Extend the Array class
Array.prototype.max = function () {
  return Math.max.apply(null, this)
}
Array.prototype.min = function () {
  return Math.min.apply(null, this)
}
Array.prototype.mean = function () {
  let i, sum
  for (i = 0, sum = 0; i < this.length; i++) sum += this[i]
  return sum / this.length
}
Array.prototype.pip = function (x, y) {
  let i
  let j
  let c = false
  for (i = 0, j = this.length - 1; i < this.length; j = i++) {
    if (
      (this[i][1] > y) !== (this[j][1] > y) &&
      x <
        ((this[j][0] - this[i][0]) * (y - this[i][1])) /
          (this[j][1] - this[i][1]) +
          this[i][0]
    ) {
      c = !c
    }
  }
  return c
}

interface IGrid {
  list: Array<Array<number>>;
  xlim: Array<number>;
  ylim: Array<number>;
  zlim: Array<number>;
  width: number;
}

interface IColor {
  min: number;
  max: number;
  color: string;
}

class VariogramClass {
  t: number[];
  x: number[];
  y: number[];
  nugget: number;
  range: number;
  sill: number;
  A: number;
  n: number;
  K: number[];
  M: number[];
  model: string;

  constructor (
    t: Array<number>,
    x: Array<number>,
    y: Array<number>,
    model: string
  ) {
    this.t = t
    this.x = x
    this.y = y
    this.nugget = 0.0
    this.range = 0.0
    this.sill = 0.0
    this.A = 1 / 3
    this.n = 0
    this.K = [] as number[]
    this.M = [] as number[]
    this.model = model
  }

  // Variogram models(变差函数模型)
  krigingVariogramGaussian (
    h: number,
    nugget: number,
    range: number,
    sill: number,
    A: number
  ) {
    return (
      nugget +
      ((sill - nugget) / range) *
        (1.0 - Math.exp(-(1.0 / A) * Math.pow(h / range, 2)))
    )
  }

  krigingVariogramExponential (
    h: number,
    nugget: number,
    range: number,
    sill: number,
    A: number
  ) {
    return (
      nugget +
      ((sill - nugget) / range) * (1.0 - Math.exp(-(1.0 / A) * (h / range)))
    )
  }

  krigingVariogramSpherical (
    h: number,
    nugget: number,
    range: number,
    sill: number
  ) {
    if (h > range) return nugget + (sill - nugget) / range
    return (
      nugget +
      ((sill - nugget) / range) *
        (1.5 * (h / range) - 0.5 * Math.pow(h / range, 3))
    )
  }

  judgeType (h: number, nugget: number, range: number, sill: number, A: number) {
    switch (this.model) {
      case 'gaussian':
        return this.krigingVariogramGaussian(h, nugget, range, sill, A)
      case 'exponential':
        return this.krigingVariogramExponential(h, nugget, range, sill, A)
      case 'spherical':
        return this.krigingVariogramExponential(h, nugget, range, sill, A)
      default:
        return 0
    }
  }
}

export default class KrigingClass {
  createArrayWithValues = (value: number, n: number) => {
    const array = []
    for (let i = 0; i < n; i++) {
      array.push(value)
    }
    return array
  };

  // Matrix algebra(矩阵代数)
  krigingMatrixDiag = (c: number, n: number) => {
    const Z = this.createArrayWithValues(0, n * n)
    for (let i = 0; i < n; i++) Z[i * n + i] = c
    return Z
  };

  krigingMatrixTranspose = (X: number[], n: number, m: number) => {
    let i
    let j
    const Z = Array(m * n)
    for (i = 0; i < n; i++) for (j = 0; j < m; j++) Z[j * n + i] = X[i * m + j]
    return Z
  };

  krigingMatrixScale = (X: number[], c: number, n: number, m: number) => {
    let i, j
    for (i = 0; i < n; i++) for (j = 0; j < m; j++) X[i * m + j] *= c
  };

  krigingMatrixAdd = (X: number[], Y: number[], n: number, m: number) => {
    let i
    let j
    const Z = Array(n * m)
    for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) Z[i * m + j] = X[i * m + j] + Y[i * m + j]
    }
    return Z
  };

  // Naive matrix multiplication(朴素矩阵乘法)
  krigingMatrixMultiply = (
    X: number[],
    Y: number[],
    n: number,
    m: number,
    p: number
  ) => {
    let i
    let j
    let k
    const Z = Array(n * p)
    for (i = 0; i < n; i++) {
      for (j = 0; j < p; j++) {
        Z[i * p + j] = 0
        for (k = 0; k < m; k++) Z[i * p + j] += X[i * m + k] * Y[k * p + j]
      }
    }
    return Z
  };

  // Cholesky decomposition(乔里斯基分解)
  krigingMatrixChol = (X: number[], n: number) => {
    let i
    let j
    let k
    const p = Array(n)
    for (i = 0; i < n; i++) p[i] = X[i * n + i]
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) p[i] -= X[i * n + j] * X[i * n + j]
      if (p[i] <= 0) return false
      p[i] = Math.sqrt(p[i])
      for (j = i + 1; j < n; j++) {
        for (k = 0; k < i; k++) X[j * n + i] -= X[j * n + k] * X[i * n + k]
        X[j * n + i] /= p[i]
      }
    }
    for (i = 0; i < n; i++) X[i * n + i] = p[i]
    return true
  };

  // Inversion of cholesky decomposition(反演乔里斯基分解)
  krigingMatrixChol2invl = (X: number[], n: number) => {
    let i, j, k, sum
    for (i = 0; i < n; i++) {
      X[i * n + i] = 1 / X[i * n + i]
      for (j = i + 1; j < n; j++) {
        sum = 0
        for (k = i; k < j; k++) sum -= X[j * n + k] * X[k * n + i]
        X[j * n + i] = sum / X[j * n + j]
      }
    }
    for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) X[i * n + j] = 0
    for (i = 0; i < n; i++) {
      X[i * n + i] *= X[i * n + i]
      for (k = i + 1; k < n; k++) X[i * n + i] += X[k * n + i] * X[k * n + i]
      for (j = i + 1; j < n; j++) {
        for (k = j; k < n; k++) X[i * n + j] += X[k * n + i] * X[k * n + j]
      }
    }
    for (i = 0; i < n; i++) for (j = 0; j < i; j++) X[i * n + j] = X[j * n + i]
  };

  // Inversion via gauss-jordan elimination(反演高斯－若尔当消元法)
  krigingMatrixSolve = (X: number[], n: number) => {
    const m = n
    const b = Array(n * n)
    const indxc = Array(n)
    const indxr = Array(n)
    const ipiv = Array(n)
    let i = 0
    let icol = 0
    let irow = 0
    let j = 0
    let k = 0
    let l = 0
    let ll = 0
    let big, dum, pivinv, temp

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (i === j) b[i * n + j] = 1
        else b[i * n + j] = 0
      }
    }
    for (j = 0; j < n; j++) ipiv[j] = 0
    for (i = 0; i < n; i++) {
      big = 0
      for (j = 0; j < n; j++) {
        if (ipiv[j] !== 1) {
          for (k = 0; k < n; k++) {
            if (ipiv[k] === 0) {
              if (Math.abs(X[j * n + k]) >= big) {
                big = Math.abs(X[j * n + k])
                irow = j
                icol = k
              }
            }
          }
        }
      }
      ++ipiv[icol]

      if (irow !== icol) {
        for (l = 0; l < n; l++) {
          temp = X[irow * n + l]
          X[irow * n + l] = X[icol * n + l]
          X[icol * n + l] = temp
        }
        for (l = 0; l < m; l++) {
          temp = b[irow * n + l]
          b[irow * n + l] = b[icol * n + l]
          b[icol * n + l] = temp
        }
      }
      indxr[i] = irow
      indxc[i] = icol

      if (X[icol * n + icol] === 0) return false // Singular

      pivinv = 1 / X[icol * n + icol]
      X[icol * n + icol] = 1
      for (l = 0; l < n; l++) X[icol * n + l] *= pivinv
      for (l = 0; l < m; l++) b[icol * n + l] *= pivinv

      for (ll = 0; ll < n; ll++) {
        if (ll !== icol) {
          dum = X[ll * n + icol]
          X[ll * n + icol] = 0
          for (l = 0; l < n; l++) X[ll * n + l] -= X[icol * n + l] * dum
          for (l = 0; l < m; l++) b[ll * n + l] -= b[icol * n + l] * dum
        }
      }
    }
    for (l = n - 1; l >= 0; l--) {
      if (indxr[l] !== indxc[l]) {
        for (k = 0; k < n; k++) {
          temp = X[k * n + indxr[l]]
          X[k * n + indxr[l]] = X[k * n + indxc[l]]
          X[k * n + indxc[l]] = temp
        }
      }
    }

    return true
  };

  // Train using gaussian processes with bayesian priors(使用高斯过程和贝叶斯先验进行训练)
  train = (
    t: Array<number>,
    x: Array<number>,
    y: Array<number>,
    model: string,
    sigma2: number,
    alpha: number
  ) => {
    const variogram = new VariogramClass(t, x, y, model)

    // Lag distance/semivariance
    let i
    let j
    let k
    let l
    let n = t.length
    const distance = Array((n * n - n) / 2)
    for (i = 0, k = 0; i < n; i++) {
      for (j = 0; j < i; j++, k++) {
        distance[k] = Array(2)
        distance[k][0] = Math.pow(
          Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2),
          0.5
        )
        distance[k][1] = Math.abs(t[i] - t[j])
      }
    }
    distance.sort(function (a, b) {
      return a[0] - b[0]
    })
    variogram.range = distance[(n * n - n) / 2 - 1][0]

    // Bin lag distance
    const lags = (n * n - n) / 2 > 30 ? 30 : (n * n - n) / 2
    const tolerance = variogram.range / lags
    const lag = this.createArrayWithValues(0, lags)
    const semi = this.createArrayWithValues(0, lags)
    if (lags < 30) {
      for (l = 0; l < lags; l++) {
        lag[l] = distance[l][0]
        semi[l] = distance[l][1]
      }
    } else {
      for (
        i = 0, j = 0, k = 0, l = 0;
        i < lags && j < (n * n - n) / 2;
        i++, k = 0
      ) {
        while (distance[j][0] <= (i + 1) * tolerance) {
          lag[l] += distance[j][0]
          semi[l] += distance[j][1]
          j++
          k++
          if (j >= (n * n - n) / 2) break
        }
        if (k > 0) {
          lag[l] /= k
          semi[l] /= k
          l++
        }
      }
      if (l < 2) return variogram // Error: Not enough points
    }

    // Feature transformation
    n = l
    variogram.range = lag[n - 1] - lag[0]
    const X = this.createArrayWithValues(1, 2 * n)
    const Y = Array(n)
    const A = variogram.A
    for (i = 0; i < n; i++) {
      switch (model) {
        case 'gaussian':
          X[i * 2 + 1] =
            1.0 - Math.exp(-(1.0 / A) * Math.pow(lag[i] / variogram.range, 2))
          break
        case 'exponential':
          X[i * 2 + 1] =
            1.0 - Math.exp((-(1.0 / A) * lag[i]) / variogram.range)
          break
        case 'spherical':
          X[i * 2 + 1] =
            1.5 * (lag[i] / variogram.range) -
            0.5 * Math.pow(lag[i] / variogram.range, 3)
          break
      }
      Y[i] = semi[i]
    }

    // Least squares
    const Xt = this.krigingMatrixTranspose(X, n, 2)
    let Z = this.krigingMatrixMultiply(Xt, X, 2, n, 2)
    Z = this.krigingMatrixAdd(Z, this.krigingMatrixDiag(1 / alpha, 2), 2, 2)
    const cloneZ = Z.slice(0)
    if (this.krigingMatrixChol(Z, 2)) this.krigingMatrixChol2invl(Z, 2)
    else {
      this.krigingMatrixSolve(cloneZ, 2)
      Z = cloneZ
    }
    const W = this.krigingMatrixMultiply(
      this.krigingMatrixMultiply(Z, Xt, 2, 2, n),
      Y,
      2,
      n,
      1
    )

    // Variogram parameters
    variogram.nugget = W[0]
    variogram.sill = W[1] * variogram.range + variogram.nugget
    variogram.n = x.length

    // Gram matrix with prior
    n = x.length
    let K = Array(n * n)
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        K[i * n + j] = variogram.judgeType(
          Math.pow(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2), 0.5),
          variogram.nugget,
          variogram.range,
          variogram.sill,
          variogram.A
        )
        K[j * n + i] = K[i * n + j]
      }
      K[i * n + i] = variogram.judgeType(
        0,
        variogram.nugget,
        variogram.range,
        variogram.sill,
        variogram.A
      )
    }

    // Inverse penalized Gram matrix projected to target vector
    let C = this.krigingMatrixAdd(K, this.krigingMatrixDiag(sigma2, n), n, n)
    const cloneC = C.slice(0)
    if (this.krigingMatrixChol(C, n)) this.krigingMatrixChol2invl(C, n)
    else {
      this.krigingMatrixSolve(cloneC, n)
      C = cloneC
    }

    // Copy unprojected inverted matrix as K
    K = C.slice(0)
    const M = this.krigingMatrixMultiply(C, t, n, n, 1)
    variogram.K = K
    variogram.M = M

    return variogram
  };

  // Model prediction(模型预测)
  predict = (x: number, y: number, variogram: VariogramClass) => {
    let i
    const k = Array(variogram.n)
    for (i = 0; i < variogram.n; i++) {
      k[i] = variogram.judgeType(
        Math.pow(
          Math.pow(x - variogram.x[i], 2) + Math.pow(y - variogram.y[i], 2),
          0.5
        ),
        variogram.nugget,
        variogram.range,
        variogram.sill,
        variogram.A
      )
    }
    return this.krigingMatrixMultiply(k, variogram.M, 1, variogram.n, 1)[0]
  };

  variance = (x: number, y: number, variogram: VariogramClass) => {
    let i
    const k = Array(variogram.n)
    for (i = 0; i < variogram.n; i++) {
      k[i] = variogram.judgeType(
        Math.pow(
          Math.pow(x - variogram.x[i], 2) + Math.pow(y - variogram.y[i], 2),
          0.5
        ),
        variogram.nugget,
        variogram.range,
        variogram.sill,
        variogram.A
      )
    }
    return (
      variogram.judgeType(
        0,
        variogram.nugget,
        variogram.range,
        variogram.sill,
        variogram.A
      ) +
      this.krigingMatrixMultiply(
        this.krigingMatrixMultiply(k, variogram.K, 1, variogram.n, variogram.n),
        k,
        1,
        variogram.n,
        1
      )[0]
    )
  };

  // Gridded matrices or contour paths(网格矩阵或轮廓路径)
  grid = (
    polygons: Array<Array<Array<number>>>,
    variogram: VariogramClass,
    width: number
  ) => {
    let i
    let j
    let k
    const n = polygons.length
    if (n === 0) return

    // Boundaries of polygons space
    const xlim = [polygons[0][0][0], polygons[0][0][0]]
    const ylim = [polygons[0][0][1], polygons[0][0][1]]
    for (
      i = 0;
      i < n;
      i++ // Polygons
    ) {
      for (j = 0; j < polygons[i].length; j++) {
        // Vertices
        if (polygons[i][j][0] < xlim[0]) xlim[0] = polygons[i][j][0]
        if (polygons[i][j][0] > xlim[1]) xlim[1] = polygons[i][j][0]
        if (polygons[i][j][1] < ylim[0]) ylim[0] = polygons[i][j][1]
        if (polygons[i][j][1] > ylim[1]) ylim[1] = polygons[i][j][1]
      }
    }

    // Alloc for O(n^2) space
    let xtarget, ytarget
    const a = Array(2)
    const b = Array(2)
    const lxlim = Array(2) // Local dimensions
    const lylim = Array(2) // Local dimensions
    const x = Math.ceil((xlim[1] - xlim[0]) / width)
    const y = Math.ceil((ylim[1] - ylim[0]) / width)

    const A: IGrid = {
      list: Array(x + 1),
      xlim: Array<number>(2),
      ylim: Array<number>(2),
      zlim: Array<number>(2),
      width: 0
    }
    for (i = 0; i <= x; i++) A.list[i] = Array(y + 1)
    for (i = 0; i < n; i++) {
      // Range for polygons[i]
      lxlim[0] = polygons[i][0][0]
      lxlim[1] = lxlim[0]
      lylim[0] = polygons[i][0][1]
      lylim[1] = lylim[0]
      for (j = 1; j < polygons[i].length; j++) {
        // Vertices
        if (polygons[i][j][0] < lxlim[0]) lxlim[0] = polygons[i][j][0]
        if (polygons[i][j][0] > lxlim[1]) lxlim[1] = polygons[i][j][0]
        if (polygons[i][j][1] < lylim[0]) lylim[0] = polygons[i][j][1]
        if (polygons[i][j][1] > lylim[1]) lylim[1] = polygons[i][j][1]
      }

      // Loop through polygon subspace
      a[0] = Math.floor(
        (lxlim[0] - ((lxlim[0] - xlim[0]) % width) - xlim[0]) / width
      )
      a[1] = Math.ceil(
        (lxlim[1] - ((lxlim[1] - xlim[1]) % width) - xlim[0]) / width
      )
      b[0] = Math.floor(
        (lylim[0] - ((lylim[0] - ylim[0]) % width) - ylim[0]) / width
      )
      b[1] = Math.ceil(
        (lylim[1] - ((lylim[1] - ylim[1]) % width) - ylim[0]) / width
      )
      for (j = a[0]; j <= a[1]; j++) {
        for (k = b[0]; k <= b[1]; k++) {
          xtarget = xlim[0] + j * width
          ytarget = ylim[0] + k * width
          if (polygons[i].pip(xtarget, ytarget)) {
            A.list[j][k] = this.predict(xtarget, ytarget, variogram)
          }
        }
      }
    }
    A.xlim = xlim
    A.ylim = ylim
    A.zlim = [variogram.t.min(), variogram.t.max()]
    A.width = width
    return A
  };
  // contour =  (value, polygons, variogram) => {};

  // Plotting on the DOM(在DOM上绘图)
  plot = (
    canvas: HTMLCanvasElement,
    grid: IGrid,
    xlim: Array<number>,
    ylim: Array<number>,
    colors: Array<IColor>
  ) => {
    // Clear screen
    const ctx = canvas.getContext('2d')
    ctx && ctx.clearRect(0, 0, canvas.width, canvas.height)
    // Starting boundaries
    const range = [
      xlim[1] - xlim[0],
      ylim[1] - ylim[0],
      grid.zlim[1] - grid.zlim[0]
    ]
    let i, j, x, y, z
    const n = grid.list.length
    const m = grid.list[0].length
    const wx = Math.ceil((grid.width * canvas.width) / (xlim[1] - xlim[0]))
    const wy = Math.ceil((grid.width * canvas.height) / (ylim[1] - ylim[0]))
    for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
        if (grid.list[i][j] === undefined) continue
        x =
          (canvas.width * (i * grid.width + grid.xlim[0] - xlim[0])) / range[0]
        y =
          canvas.height *
          (1 - (j * grid.width + grid.ylim[0] - ylim[0]) / range[1])
        z = (grid.list[i][j] - grid.zlim[0]) / range[2]
        if (z < 0.0) z = 0.0
        if (z > 1.0) z = 1.0
        if (ctx) {
          ctx.fillStyle = this.getColor(colors, grid.list[i][j])
        }
        ctx &&
          ctx.fillRect(Math.round(x - wx / 2), Math.round(y - wy / 2), wx, wy)
      }
    }
  };

  // custom color(自定义色彩)
  getColor = (colors: Array<IColor>, z: number) => {
    const l = colors.length
    for (let i = 0; i < l; i++) {
      if (z >= colors[i].min && z < colors[i].max) return colors[i].color
    }
    if (z < 0) {
      return colors[0].color
    } else {
      return ''
    }
  };
}
