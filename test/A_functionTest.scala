package UnbiasedEstimator

import org.junit._
import Assert._

import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.noise.ComplexGaussian
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util
import numbers.finite.integration.RealIntegral.trapezoidal

class A_functionTest {
   
  /**  @Before
  def setUp: Unit = {
  }  **/

 /** @After
  def tearDown: Unit = {
  }  **/

  @Test
  def testA_function = {
    
  final def sin(x : Double) = scala.math.sin(x)
  final def cos(x : Double) = scala.math.cos(x) 
  final def abs(x : Double) = math.abs(x)
  final def sqr(x : Double) = x*x
  final def cub(x : Double) = x*x*x
  final def sqrt(x : Double) = scala.math.sqrt(x) 
  final def exp(x : Double) = scala.math.exp(x)
  final def erf(x : Double) = pubsim.Util.erf(x)
  val pi = scala.math.Pi  
  
  val M = 4 
  val L = 32 
  val p = 0.5
  val d = 0.5    
  
   def Bisection(A: Double => Double, rhohat: Double, a: Double, b : Double, tolerance: Double): Double = { 
    var aa = a
    var bb = b
    while (Math.abs(A(bb)-rhohat) > tolerance){
      var c = (aa + bb)/2.0 
      if ( (A(c)<rhohat) )
        aa = c
      else
        bb = c
    } 
  return bb    
    } 
    
val SNRdBs = -20 to 20 by 1
val SNRs = SNRdBs.map(db => scala.math.pow(10.0, db/10.0))
val noises = SNRs.map( snr => new ComplexGaussian(a0.mag2/snr/2.0) ) 
val sigma = 0.01  
 

    
  def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi) 

  def g(phi : Double, rho0 : Double, sigma : Double) : Double = { 
  val a = cos(phi)  
  val b = sin(phi)
  val k =  rho0/sigma // rho0/sigma
  val PHI = (1 + erf(a*k/sqrt(2)))/2
  
  return a*exp(-k*k/2)/2/pi + 1/k/sqrt(2*pi) * exp(-k*k/2*b*b) * PHI * (1 + k*k*a*a) 
  }
    def h1( rho0:Double, sigma:Double) = trapezoidal ( phi => cos(0 + phi)*g(phi, rho0 : Double, sigma : Double), -pi, pi, 200 )
    def h2( rho0:Double, sigma:Double) = trapezoidal ( phi => cos(fracpart(0 + phi))*g(phi,rho0 : Double, sigma : Double), -pi, pi, 200 )
    def G(  rho0:Double, sigma:Double) = p*h1(rho0 : Double, sigma : Double) + d*h2(rho0 : Double, sigma : Double)
     
    def A(rho0:Double)= trapezoidal( phi => rho0*G( rho0, sigma ), -pi, pi, 200)     
   
  
  println(A(1))
  
  
  }
  
    
    
    
}
