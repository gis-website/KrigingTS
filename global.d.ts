/*
 * @Author: TQtong 2733707740@qq.com
 * @Date: 2023-05-08 20:21:56
 * @LastEditors: TQtong 2733707740@qq.com
 * @LastEditTime: 2023-07-27 14:08:35
 * @FilePath: \KrigingTS\global.d.ts
 * @Description: global
 */

declare global {
    interface Array {
        max():number,
        min():number,
        mean():number,
        pip(x:number, y:number):boolean
    }
}
