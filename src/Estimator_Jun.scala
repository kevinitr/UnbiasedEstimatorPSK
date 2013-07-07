package Estimator

import cam.psk.MackenthunCoherent
import cam.psk.ComplexAmplitudeEstimator
import cam.noise.ComplexGaussian
import pubsim.distributions.complex.SymmetricComplexNormal
import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import pubsim.distributions.complex.ComplexRandomVariable
import pubsim.Util
import numbers.finite.integration.RealIntegral.trapezoidal
 
object UnbiasedEstimator  {
  def main(args: Array[String]) {}      
  final def sin(x : Double) = scala.math.sin(x)
  final def cos(x : Double) = scala.math.cos(x) 
  final def abs(x : Double) = math.abs(x)
  final def sqr(x : Double) = x*x
  final def cub(x : Double) = x*x*x
  final def sqrt(x : Double) = scala.math.sqrt(x) 
  final def exp(x : Double) = scala.math.exp(x)
  final def erf(x : Double) = pubsim.Util.erf(x)
  val pi = scala.math.Pi  
 
   val Ms = List(2, 4) 
   val Ls = List( 256, 1024, 4096 )    
 
 def Bisection(A: Double => Double, rhohat: Double, a: Double, b : Double, tolerance: Double): Double = { 
    var aa = a
    var bb = b
    while (abs(A(bb)-rhohat) > tolerance){
      var c = (aa + bb)/2.0 
      if ( (A(c)<rhohat) )
        aa = c
      else
        bb = c
    } 
 return bb    
    } 
val a0 = new PolarComplex(1,2*scala.math.Pi*(new scala.util.Random).nextDouble) //set rho0=1
val iters = 1000  // 1000 replications
val m = 2
 
val SNRdBs = -10 to 10 by 1
val SNRs = SNRdBs.map( db => scala.math.pow(10.0, db/10.0) )
val noises = SNRs.map( snr => new ComplexGaussian(a0.mag2/snr/2.0) )
val sigma = SNRs.map( snr => sqrt(a0.mag2/snr/2.0) )    
 
for( L <- Ls; M <- Ms ) {
  for( numpilots <- List( L/16, L/2  ).distinct.filterNot(b => b == 0) ) {
    val P = 0 until numpilots //pilots at the front
    val D = numpilots until L //data at the back
    val p = P.length.toDouble/L 
    val d= 1-p 
 
  def fracpart(x : Double) = x - 2*pi/M*scala.math.round(M*x/2/pi) 

  def g(phi : Double, rho0 : Double, sigma : Double) : Double = { 
  val a = cos(phi)  
  val b = sin(phi)
  val k =  rho0/sigma  
  val PHI = (1 + erf(a*k/sqrt(2)))/2 
  
  return a*exp(-k*k/2)/2/pi + 1/k/sqrt(2*pi) * exp(-k*k/2*b*b) * PHI * (1 + k*k*a*a) 
  }
    def h1( rho0:Double, sigma:Double) = trapezoidal ( phi => cos(0 + phi)*g(phi, rho0 : Double, sigma : Double), -pi, pi, 200 )
    def h2( rho0:Double, sigma:Double) = trapezoidal ( phi => cos(fracpart(0 + phi))*g(phi,rho0 : Double, sigma : Double), -pi, pi, 200 )
    def G(  rho0:Double, sigma:Double) = p*h1(rho0 : Double, sigma : Double) + d*h2(rho0 : Double, sigma : Double) 
 
    def estfactory = List( 
      (p : Seq[Complex]) => new MackenthunCoherent(M,P,D,p)
    )
    for( estf <- estfactory ) {    
   //  val estname =  estf(null).getClass.getSimpleName + "M" + M + "L" + L.toString + "absP" + P.length
     val estname =    "M" + M + "L" + L.toString + "P" + P.length
     print("Running " + estname +"\n")      
     val mselist = noises.indices.map { i =>	   
  
	val rand = new scala.util.Random
	val s = (1 to L).map(m => new PolarComplex(1, 2*scala.math.Pi*rand.nextInt(M)/M)) 
	val est = estf(s)  	  
        val p = P.length.toDouble/L
        var msebiased  = 0.0; var mseunbiased  = 0.0;  var mseunbiasedsigma=0.0
        var msesigma = 0.0; var sigma2 = 0.0; 
    
    // Estimate variance of noise from a period of silence       
        for( itr <- 1 to m ) {
    val pilot = 1.0 to L/m by 1.0     // use same length for each subinterval as received data block 
    val pilotlength = pilot.length
    val noisesreal = pilot.map( p => noises(i).noise.real )   // generate noise 
    val noisesimag = pilot.map( p => noises(i).noise.imag ) 
    val variancereal = noisesreal.foldLeft(0.0)((r,v) => ( r + sqr(v) )) /pilotlength
    val varianceimag = noisesimag.foldLeft(0.0)((r,v) => ( r + sqr(v) )) /pilotlength
    sigma2 += sqrt((variancereal+varianceimag)/2) 
    }          
    
     for( itr <- 1 to iters ) {
 
     val y = s.map( si => a0*si + noises(i).noise)
     val ahat = est.estimate(y)  
 
     def A(rho0:Double)=  rho0*G( rho0 , sigma(i))    
     def A1(rho0:Double)=  rho0*G( rho0 , sigma2 / m )  //Simulation with estimated sigma in real time
     msebiased += (ahat - a0).mag2    
    var rho1= Bisection(A,ahat.magnitude,0.0, 10.0, 0.001) 
    var rho2= Bisection(A1,ahat.magnitude,0.0, 10.0, 0.001)
     mseunbiased += ( rho1 - a0.magnitude )*( rho1 - a0.magnitude )  
     mseunbiasedsigma += ( rho2 - a0.magnitude )*( rho2 - a0.magnitude )  
   }  
            
    print(msebiased /iters +"\t" +  mseunbiased /iters +"\t" + mseunbiasedsigma/iters+ "\n")
    (msebiased / iters, mseunbiased / iters, mseunbiasedsigma / iters)   
      }.toList 
      
         val filea = new java.io.FileWriter("data/" + estname + "a")
         val fileb = new java.io.FileWriter("data/" + estname + "b")
         val filec = new java.io.FileWriter("data/" + estname + "c")
          (mselist, SNRdBs).zipped.foreach{ (mse, snr) =>
        val (msebiased, mseunbiased,mseunbiasedsigma) = mse
 
         filea.write(snr.toString.replace('E', 'e') + "\t" +  msebiased  + "\n") 
         fileb.write(snr.toString.replace('E', 'e') + "\t" +   mseunbiased  + "\n") 
         filec.write(snr.toString.replace('E', 'e') + "\t" +   mseunbiasedsigma  + "\n") 
     }
        filea.close;fileb.close;filec.close}  

     }           
   }   
 }   