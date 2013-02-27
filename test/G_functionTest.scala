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

import org.junit._
import Assert._
 
class G_functionTest {
  final def sin(x : Double) = scala.math.sin(x)
  final def cos(x : Double) = scala.math.cos(x) 
  final def abs(x : Double) = math.abs(x)
  final def sqr(x : Double) = x*x
  final def cub(x : Double) = x*x*x
  final def sqrt(x : Double) = scala.math.sqrt(x) 
  final def exp(x : Double) = scala.math.exp(x)
  final def erf(x : Double) = pubsim.Util.erf(x)
  
  /**  @Before
  def setUp: Unit = {
  }  **/

 /** @After
  def tearDown: Unit = {
  }  **/

  @Test def G_function={  //test G function by chaging SNR
    
 
  val pi = scala.math.Pi  
  val M = 2 
  val Ls = 32
  val p = 0.5 
  val d = 0.5

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
    
  var G1=G(1,1)
  var G2=G(10,10)
  
  assertEquals(G1,G2)
  println("G_functionTest ");
    
  }
}
