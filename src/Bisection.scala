/* Bisection method */
package UnbiasedEstimator

class Bisection {



  def Bisection(A: Double => Double, rhohat: Double, a:Double, b:Double, tolerance:Double): Double = {
    var aa = a
    var bb = b
    while (Math.abs(A(bb)-rhohat) > tolerance)
    {
      var c = (aa + bb)/2.0
      if ( (A(c)<rhohat))
        bb = c
      else
        aa = c
    }
    return aa
  }
}
