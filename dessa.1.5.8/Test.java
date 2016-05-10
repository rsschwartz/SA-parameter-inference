
import java.util.Random;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

/**
 * The <code>Test</code> class loads a <code>Simulation</code> 
 * from XML, starts it and lists important control variables. 
 * 
 * @author Rori Rohlfs
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @version 1.4
 */
public class Test  {

	/* ******************** Output Modifiers ********************/

	/** Controls the number of <code>Event</code>s between outputs */
	public static int eventsPerPrint = 1000;

	/** Weather or not to print the <code>Conformation</code> distribution */
	public static boolean printConfDist = false;

	/** Maximum size for which distribution statistics are printed */
	public static int maxOutputSize = 8;

	/* ******************** General Values ********************/

	/** The name for this Simulation type */
	public static String simType = "xml";

	/** Which xml file to use */
	public static String mode = "";

	/** Mass for a <code>Subunit</code> */
	public static double subMass = 3.4;

	/** Ratio of <code>Subunit</code> to <code>BindingSite</code>, should be <=0.5 */
	public static double subunitToBS = 0.5;

	/** Height of a <code>BindingSite</code> */
	public static float bindingSiteHeight = 0.10f;

	/** Radius for a <code>Subunit</code> */
	public static double subRadius = bindingSiteHeight * subunitToBS;

	/* ******************** Simulation Variables ********************/

	/** Use the diffusion model */
	//public static boolean diffusionEnabled = false;

	/** Pick events that only allow monomer binding */
	public static boolean bindMonomerOnly = false;

	/** only sample breaking event for bonds not involving a loop */
	public static boolean noLoopOnly = true;

	/** For filament assembly, allow only monomers break off */
	public static boolean breakOnlyEnds = false;

	/** Enable conformational switching event sampling */
	public static boolean csAllowed = true;

	/** The int returned from Assemblies.numbSubunits() that will be treated as a monomer */
	public static int sizeOfSubunit = 1;

	/**  Maximum allowed assembly size */
	public static int maxLength = 360;

	/** The amount of volume to search for nearest neighbors */
	public static double binSize = 0.1;

	/** The distance tolerance for loop detection*/
	public static double distanceTolerance = 0.1;

	/** Controls weather or not the Spring Force model is used */
	public static boolean springForce = false;

	/* ******************** Random Number Generator ********************/
	
	/** The Random number generator for the simulation */
	public static Random rand;

	/** The random number generator for the seed */
	public static Random seedGenerator;

	/** The seed used for the random number generator */
	public static int seedUsed;

	/** MAx Simulation time allowed. After which simulation turns off automatically */
	public static double maxSimulationTime;

	/**
	 * main method to start the simulator
	 */
	public static void main(String args[]) {

		if (args.length == 0) {
			System.out.println("usage: xmlFile.xml [Events Per Printout]" +
			"[Max Simulation Time] [Random Seed] [Maximum Output Size] ");
			System.exit(-1);
		}

		mode = args[0];
		
		try {
			eventsPerPrint = Integer.parseInt(args[1]);
		} catch (Exception e) {
			eventsPerPrint = 10000;
		}

		try {
			maxSimulationTime = Double.parseDouble(args[2]);
		} catch (Exception e ) {
			maxSimulationTime =  Double.MAX_VALUE;
		}
		try {
			seedUsed = Integer.parseInt(args[3]);
		} catch (Exception e ) {
			seedGenerator = new Random();
			seedUsed = seedGenerator.nextInt();
		}
		
		try {
			maxOutputSize = Integer.parseInt(args[4]);
		} catch(Exception e) {
			maxOutputSize = 200;
		}
		
		
		System.out.println(String.format("SeedUsed: %d", seedUsed));
		rand = new Random(seedUsed);     


		System.out.println("Self-assembly Simulation (Dessa v1.5.8)");

		XMLReader reader = null;
		try {
			reader = new XMLReader(mode);
		} catch (Exception e) {
			e.printStackTrace();
		}

		Simulation sim  = reader.getSim();
		sim.run();

	}

	/**
	 * Rotates a Vector around an axis by the specified angle.
	 * 
	 * @param v - The vector to rotate
	 * @param a - the AxisAngle to rotate around
	 * @return Vector3d - the rotated Vector3d
	 */
	public static Vector3d rotateByAxis(Vector3d v, AxisAngle4d a) {

		Matrix4d m = new Matrix4d();
		m.set(a);
		Vector3d rv = new Vector3d();
		m.transform(v, rv);
		return rv;
	}

	/**
	 * 
	 * @param axis
	 * @param upV
	 * @return
	 */
	public static Vector3d makePerpendicular(Vector3d axis, Vector3d upV) {

		double angle = axis.angle(upV);

		Vector3d projectedV = new Vector3d();

		if (angle < 0.00001 || (Math.PI-angle)<0.00001 ) {
			System.out.println("Bad_choice_in_upVector");
			System.exit(0);
		} else {

			Vector3d tem = new Vector3d(axis);
			tem.scale(axis.dot(upV) / axis.lengthSquared());

			projectedV.sub(upV, tem);
		}
		projectedV.normalize();

		return projectedV;
	}
}