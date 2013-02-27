package UnbiasedEstimator

import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.noise.ComplexGaussian
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util
import numbers.finite.integration.RealIntegral.trapezoidal
 
class G_function  {
 // def main(args: Array[String]) {}      
  final def sin(x : Double) = scala.math.sin(x)
  final def cos(x : Double) = scala.math.cos(x) 
  final def abs(x : Double) = math.abs(x)
  final def sqr(x : Double) = x*x
  final def cub(x : Double) = x*x*x
  final def sqrt(x : Double) = scala.math.sqrt(x) 
  final def exp(x : Double) = scala.math.exp(x)
  final def erf(x : Double) = pubsim.Util.erf(x)
  val pi = scala.math.Pi  
  val Ms = List(2) //BPSK, QPSK, 8-PSK val Ms = List(2,4,8)
  val Ls = List(32)
 
 
for( L <- Ls; M <- Ms ) {
  for( numpilots <- List( L/2 ).distinct.filterNot(b => b == 0) ) {
    val P = 0 until numpilots //pilots at the front
    val D = numpilots until L //data at the back
    val p = P.length.toDouble/L 
    val d= 1-p 
    
  def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi) 

  def g(phi : Double, rho0 : Double, sigma : Double) : Double = { 
  val a = cos(phi)  
  val b = sin(phi)
  val k =  rho0/sigma // rho0/sigma
  val PHI = (1 + erf(a*k/sqrt(2)))/2
  
  return a*exp(-k*k/2)/2/pi + 1/k/sqrt(2*pi) * exp(-k*k/2*b*b) * PHI * (1 + k*k*a*a) 
  }
    def h1( rho0:Double, sigma:Double) = trapezoidal ( phi =>  cos(0 + phi)*g(phi, rho0 : Double, sigma : Double), -pi, pi, 200 )
    def h2(  rho0:Double, sigma:Double) = trapezoidal  ( phi => cos(fracpart(0 + phi))*g(phi,rho0 : Double, sigma : Double), -pi, pi, 200 )
    def G(  rho0:Double, sigma:Double) = p*h1(rho0 : Double, sigma : Double) + d*h2(rho0 : Double, sigma : Double)
     
      println(G(100,10))

  }
  }
  }
  
