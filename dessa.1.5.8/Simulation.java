
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;
import java.util.HashMap;
import java.util.Iterator;

import javax.vecmath.Vector3d;


//TODO Document methods and class variables
//TODO clean up code
//TODO change naming of get/sample methods to be consistent
//TODO understand the the update count methods
//TODO make it so size distribution is constantly tracked not updated at intervals

/**
 * This class runs the simulation. It picks events using a stochastic method,
 * stores them in a queue and then executes them. 
 * 
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @version 1.4
 */
public class Simulation {

	/** An array of the current size distribution in the simulation */
	private int[] currentSizeDistribution=new int[Test.maxLength];

	/** 
	 * Counter for creating new Assemblies. 
	 */
	private int assemblyNumber;


	/** Current simulation step */
	private static int currentStep = 0;

	/** Current simulation time */
	private double curtime;

	/** Solution simulation occurs in */
	//private Solution mysoln;

	/** HashMap of Assembly name to Assembly */
	private HashMap<Integer, Assembly> assembliesHashMap;

	/** The priority queue where events are stored */
	private PriorityQueue pq;

	/** 
	 * A map of the BindingTimes. Maps from name to a HashMap 
	 * of partner names and then a double[3], of bond, break, and fast bind 
	 * average times.
	 */

	private static HashMap<String, BindingSiteType> bstMap;

	/**
	 * A HashMap that maps from Conformation names to a map of
	 * potential switches of Name then Double of average switch time.
	 */
	private HashMap<String, HashMap<String, Double>> confTimes;
	private HashMap<String, Integer> freeBSs;
	private HashMap<String, Integer> monomerBSs;
	
	private HashMap<String, Integer> typemap = new HashMap<String, Integer>();//TODO Remove Later
	Integer fileCounter = new Integer(1);//TODO Remove Later
	boolean got73 = false;
	boolean got80 = false;
	/**
	 * A constructor. This is meant to be used from the XML loader.
	 * 
	 * @param initialAssemblies - Vector of all initial assemblies
	 * @param startTime - the starting time of this simulation
	 * @param bondTimes - Map of Bind/Break/FastBind times
	 * @param confTimes - Conformational Switching Map
	 * @param bindingPartner - Map of potential binding Partners
	 * @param mySolution - Solution simulation occurs in
	 */
	public Simulation(Vector<Assembly> initialAssemblies, double startTime, 
			HashMap<String, BindingSiteType> bindingSiteTypeMap,
			HashMap<String, HashMap<String, Double>> cnfTimes/*,
			Solution mySolution*/) {

		assemblyNumber = 0;

		assemblyNumber = initialAssemblies.size();
		++assemblyNumber;
		//mysoln = mySolution;

		bstMap = bindingSiteTypeMap;
		curtime = startTime;
		confTimes = cnfTimes;

		pq = new BinaryHeap();

		/** 
		 * Initializes the count HashMaps 
		 */


		freeBSs = new HashMap<String, Integer>(bstMap.size());
		monomerBSs = new HashMap<String, Integer>(bstMap.size());
		Iterator<String> it = bstMap.keySet().iterator();

		while (it.hasNext()) {
			String name = it.next();
			freeBSs.put(name, new Integer(0));
			monomerBSs.put(name, new Integer(0));
		}

		/**
		 * Stores the correct values into the count HashMaps
		 */
		int size = initialAssemblies.size();
		assembliesHashMap = new HashMap<Integer, Assembly>(size);
		for (int i = 0; i < size; ++i) 
		{
			Assembly tem = initialAssemblies.get(i);

			if(tem==null)
				System.out.println("assem is null");
			assembliesHashMap.put(tem.getID(), tem);
			HashMap<String, Integer> a = tem.getBSCounts();


			if (tem.numSubunits() == Test.sizeOfSubunit) 
				updateCounts(monomerBSs, a, true);

			updateCounts(freeBSs, a, true);

		}
		
		typemap.put("bst0a",1);
		typemap.put("bst0b",2);
		typemap.put("bst0c",3);
		typemap.put("bst0d",4);
		typemap.put("bst0e",5);
		typemap.put("bst0f",6);
		typemap.put("bst1a",7);
		typemap.put("bst1b",8);
		typemap.put("bst1c",9);
		typemap.put("bst1d",10);
		typemap.put("bst1e",11);
	
		initialize();
	}



	/**
	 * This method populates the queue with the initial events.
	 * It considers FormBond, BreakBond and if desired ConfChange events.
	 */
	private void initialize() {

		Vector<Assembly> allAssemblies = new Vector<Assembly>(assembliesHashMap.values());
		int size = allAssemblies.size();
		if (size == 0)
			return;

		Event[] formBondArray = new Event[allAssemblies.size()];

		/**
		 * Populate array with the minimum FormBondEvent for each
		 * assembly
		 */
		for (int i = 0; i < size; ++i) {

			Assembly assememblyFirst = allAssemblies.get(i);

			for (int j = 0; j < i; ++j) {

				Assembly assememblySecond = allAssemblies.get(j);

				Event tem = sampleTwoAssems(assememblyFirst,
						assememblySecond);

				double eventTime;

				if (tem == null)
					eventTime = Double.MAX_VALUE;
				else
					eventTime = tem.getEndTime();

				if (formBondArray[i] == null) 
					formBondArray[i] = tem;

				if (formBondArray[j] == null) 
					formBondArray[j] = tem;

				if (formBondArray[i] != null
						&& eventTime < formBondArray[i].getEndTime()) 
					formBondArray[i] = tem;

				if (formBondArray[j] != null
						&& eventTime < formBondArray[j].getEndTime()) 

					formBondArray[j] = tem;
			}
		}

		//breaking event and conf change event
		for (int k = 0; k < size; ++k) {

			Assembly assemi = (Assembly) allAssemblies.get(k);
			Event minFormBondEvent = formBondArray[k];

			Event brkEvt = sampleBreakBondEvent(assemi);
			Event minEvent = screenEvents(minFormBondEvent, brkEvt);

			//cnfChngEvt could be null
			if (Test.csAllowed) {
				Event cnfChngEvt = getConfChangeEvent(assemi);
				minEvent = screenEvents(cnfChngEvt, minEvent);
			}

			formBondArray[k] = minEvent;
		}
		//	remove duplicated formBondEvent
		for (int k = 0; k < size; ++k) {

			/* for virus shell assembly, there will be different bs, consider it
			 * later
			 */

			if (formBondArray[k] == null)
				continue;

			boolean repeat = false;
			for (int m = k + 1; m < size; ++m) {

				if ((formBondArray[k]) == (formBondArray[m])) {
					repeat = true;
					break;
				}
			}

			if ((repeat == false)
					&& (formBondArray[k].getEndTime() < Double.MAX_VALUE)) {				
				sendEvent(formBondArray[k]);
			}            
		}

		for(int i =0; i < size ; ++i)
			allAssemblies.get(i).setValidTime(curtime);
	}



	/**
	 * This method runs the simulation. It steps for as long as there are
	 * valid events. It will output the distribution and if desired an xml
	 * at specified events. 
	 */
	public void run() {

		int eventsPerPrintout = Test.eventsPerPrint;
		int i = 1;

		updateCount();
		printDistributionToScreen();

		while (!pq.isEmpty()) {

			step();
			i++;
			
			if(curtime >= Test.maxSimulationTime)//stop simulation after max time alloted
				break;

			if( i % eventsPerPrintout == 0 ) {
				//printToFile();
				updateCount();
				printDistributionToScreen();			
			}
		}

		updateCount();
		printDistributionToScreen();	
		//printToFile();
	}
	
	/**
	 * This keeps track of the counts of each assembly size.
	 */
	private void updateCount() 
	{
		for (int size = 0; size < currentSizeDistribution.length; ++size) 
			currentSizeDistribution[size]=0;

		Iterator<Assembly> aItr = assembliesHashMap.values().iterator();
		while(aItr.hasNext())
		{
			Assembly temA = aItr.next();
			int k = temA.numSubunits();
			currentSizeDistribution[k-1]=currentSizeDistribution[k-1]+1;
		}		
	}



	/**
	 * This method performs the actual step. It gets the next event and then
	 * sends it to be processed if it is still valid.
	 */
	public void step() {

		++currentStep;
		Event ev = (Event) pq.remove();

		curtime = ev.getEndTime();

		double[] validTimes = { -1, -1 };

		Assembly[] assemsInvolved = ev.getAssembliesInvolved();
		Assembly asm1 = assemsInvolved[0];

		//only one Assembly involved and it no longer exists
		if (assemsInvolved[1] == null && !assembliesHashMap.containsKey(asm1.getID())) 
			return;

		else if (assemsInvolved[1] != null) { //two Assemblies involved

			Assembly asm2 = (Assembly) assemsInvolved[1];

			//both Assemblies no longer exist
			if (!assembliesHashMap.containsKey(asm1.getID()) && !assembliesHashMap.containsKey(asm2.getID())) 
				return;

			//only asm1 no longer exists
			if (!assembliesHashMap.containsKey(asm1.getID()) && assembliesHashMap.containsKey(asm2.getID())) {

				//sample new events for asm2 before returning
				Assembly[] tmp = { assemsInvolved[1] };
				validTimes[1] = asm2.validTime();

				if (ev.isValid(validTimes)) 
					getNewEvents(tmp);

				return;
			}

			//only asm2 no longer exists
			if (assembliesHashMap.containsKey(asm1.getID()) && !assembliesHashMap.containsKey(asm2.getID())) {

				Assembly[] tmp = { assemsInvolved[0] };
				validTimes[0] = asm1.validTime();

				if (ev.isValid(validTimes)) 
					getNewEvents(tmp);

				return;
			}
		}

		//track validTimes
		validTimes[0] = asm1.validTime();
		if (assemsInvolved[1] != null)
			validTimes[1] = ((Assembly) assemsInvolved[1]).validTime();

		/*
		 * In the case that Assembly(s) involved in one Event exist, If event is
		 * valid, processes it, generates new events for next iteration. else
		 * picks another event for an Assembly that is valid(try again)
		 *  
		 */
		if (ev.isValid(validTimes)) {
			processEvent(ev);

			//It is tested when only at most two Assembly(s) are involved
			//One Assembly is invalid, consider the other one in case of
			// FormBondEvent
		} else {

			//Consider only maximum two Assembly(s) in any event
			if (ev.getEventType() == EventType.formBndEvt) {

				if (ev.getPostTime() >= validTimes[0]) {
					Assembly[] tmp = { assemsInvolved[0] };
					getNewEvents(tmp);

				} else if (ev.getPostTime() >= validTimes[1]) {
					Assembly[] tmp = { assemsInvolved[1] };
					getNewEvents(tmp);
				}
			}
		}
	}

	/**
	 * This method process an event. It will perform the event
	 * and pick a new one for the assemblies involved.
	 * 
	 * @param ev - Event to perform
	 */
	private void processEvent(Event ev) {

		if (ev.getEventType() == EventType.cnfChngEvt) {

			Assembly assem = ev.getAssembliesInvolved()[0];

			updateCounts(freeBSs, assem.getBSCounts(), false);
			if (assem.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(monomerBSs, assem.getBSCounts(), false);

			}

			assem.setValidTime(curtime);
			ev.getSubunit().changeConf(ev.getDomainID(), ev.getNewConf());
			ev.setAssembliesInvolved(new Assembly[] { ev.getAssembliesInvolved()[0] });

			updateCounts(freeBSs, assem.getBSCounts(), true);
			if (assem.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(monomerBSs, assem.getBSCounts(), true);
			}

		} else if (ev.getEventType() == EventType.brkBndEvt) {			
			Assembly oldasm = ev.getBS().getAssembly();
			BindingSite bs = ev.getBS();
			BindingSite partner = ev.getPartner();

			simpleUpdateCountsAll(bs, true);
			simpleUpdateCountsAll(partner, true);

			//In general, splitAssembly will produce a new Assembly            
			Assembly newasm = oldasm.splitAssembly(ev, assemblyNumber);

			//if the splitAssembly makes a new assembly
			if (newasm != null) 
			{
				assembliesHashMap.put(newasm.getID(), newasm);

				ev.setAssembliesInvolved(new Assembly[] { oldasm, newasm });
				++assemblyNumber;


				if (newasm.numSubunits() == Test.sizeOfSubunit) 
					updateCounts(monomerBSs, newasm.getBSCounts(), true);

				if (oldasm.numSubunits() == Test.sizeOfSubunit) 
					updateCounts(monomerBSs, oldasm.getBSCounts(), true);

				newasm.setValidTime(curtime);
				oldasm.setValidTime(curtime);

			} else {

				double fastBindingInterval = 0.0;
				oldasm.setValidTime(curtime);
				Assembly[] sameAssembly = { ev.getAssembliesInvolved()[0],
						ev.getAssembliesInvolved()[0] };
				Event fastBindingEvent = new Event(curtime,
						curtime + fastBindingInterval, sameAssembly, ev.getBS(),
						ev.getPartner(),EventType.formBndEvt);
				sendEvent(fastBindingEvent);

				// return to prevent picking of new events for the assembly
				return;
			}
		} else if (ev.getEventType() == EventType.formBndEvt) 
		{
			BindingSite bs = ev.getBS();
			Assembly asm1 = bs.getAssembly();

			BindingSite partner = ev.getPartner();
			Assembly asm2 = partner.getAssembly();

			if (asm1 == asm2) {

				asm1.setValidTime(curtime);
				asm1.fastBind(ev);
				ev.setAssembliesInvolved(new Assembly[] { asm1 });

				simpleUpdateCountsAll(bs, false);
				simpleUpdateCountsAll(partner, false);


			} else if ((asm1.numSubunits() + asm2.numSubunits()) <= Test.maxLength) {


				updateCounts(freeBSs, asm1.getBSCounts(), false);
				updateCounts(freeBSs, asm2.getBSCounts(), false);
				if (asm1.numSubunits() == Test.sizeOfSubunit)
					updateCounts(monomerBSs, asm1.getBSCounts(), false);

				if (asm2.numSubunits() == Test.sizeOfSubunit)
					updateCounts(monomerBSs, asm2.getBSCounts(), false);

				Assembly tmp = asm1.bindAssembly(asm2, ev);

				if (tmp != null) {
					assembliesHashMap.remove(asm2.getID());
					ev.setAssembliesInvolved(new Assembly[] { asm1 });
					asm1.setValidTime(curtime);

					updateCounts(freeBSs, asm1.getBSCounts(), true);

				} else {
					asm1.setValidTime(curtime);
					asm2.setValidTime(curtime);
					updateCounts(freeBSs, asm1.getBSCounts(), true);
					updateCounts(freeBSs, asm2.getBSCounts(), true);
					if (asm1.numSubunits() == Test.sizeOfSubunit) {
						updateCounts(monomerBSs, asm1.getBSCounts(), true);
					}
					if (asm2.numSubunits() == Test.sizeOfSubunit) {
						updateCounts(monomerBSs, asm2.getBSCounts(), true);
					}

				}
			} else { // if the assembly is too big
				asm1.setValidTime(curtime);
				asm2.setValidTime(curtime);
			}
		}

		/*
		 * generates new events for queue based on new current state, only need
		 * to update nextEvent for assemblies that were involved in "ev".if only
		 * one assembly is involved, second element is null.
		 */
		getNewEvents(ev.getAssembliesInvolved());
	}



	/**
	 * Returns the Coefficient 
	 * 
	 * @param primary
	 * @param secondary
	 * @return double - the Coefficient
	 */
	double getCoef(int primary, int secondary) {

		return Math.sqrt((primary + secondary)
				/ (2.0 * primary * secondary));
	}

	public final static Vector3d getBSTPostion(String type)
	{
		return new Vector3d(bstMap.get(type).getBSTPosition());
	}

	public final static double getBSTBindAngle(String type,String partnerType)
	{
		return bstMap.get(type).getBindingAngle(partnerType);
	}

	/**
	 * Picks new events for the two Assemblies. Puts the new Events into
	 * the queue. Events are picked that are either a FormBondEvent, 
	 * BreakBondEvent or if allowed a ConfChangeEvent. 
	 * 
	 * @param assems - Array of the two Assemblies to pick new events for
	 */
	private void getNewEvents(Assembly[] assems) {

		HashMap otherFreeBSs = null;
		HashMap<String, Integer> otherMonomerBSs = null;


		otherFreeBSs = new HashMap<String, Integer>(freeBSs);
		otherMonomerBSs = new HashMap<String, Integer>(monomerBSs);

		HashMap selfBSs0 = null;
		HashMap selfBSs1 = null;
		Event myMinEvent0 = null;
		Event myMinEvent1 = null;

		//Some error checking

		if (assems.length == 1) {

			if (assems[0] == null
					|| (this.assembliesHashMap.get(assems[0].getID()) == null)) {
				System.out
				.println("argument_error_1 in Simulation getNewEvents");
				System.exit(1);
			}

		} else if (assems.length == 2) {
			if ((assems[0] == null)
					|| (this.assembliesHashMap.get(assems[0].getID()) == null)
					|| (assems[1] == null)
					|| (this.assembliesHashMap.get(assems[1].getID()) == null)) {
				System.out
				.println("argument_error_2 in Simulation getNewEvents");
				System.exit(1);
			}
			if (assems[0].equals(assems[1])) {
				System.out
				.println("argument_error_3 in Simulation getNewEvents");
				System.exit(1);
			}
		} else {
			System.out.println("argument_error_4 in Simulation getNewEvents");
			System.exit(1);
		}

		if (assems.length == 1) {

			Assembly a0 = (Assembly) assems[0];
			selfBSs0 = a0.getBSCounts();

			updateCounts(otherFreeBSs, selfBSs0, false);
			a0.setValidTime(curtime);
			myMinEvent0 = getFormBondEvent(a0, otherFreeBSs,
					otherMonomerBSs, null);

			myMinEvent0 = screenEvents(myMinEvent0, getConfChangeEvent(a0));
			Event brkEvt = sampleBreakBondEvent(a0);
			myMinEvent0 = screenEvents(myMinEvent0, brkEvt);


			if (a0.numSubunits() == Test.sizeOfSubunit) {
				if (brkEvt != null)
					System.out.println("error1insimulation");
			}


		} else {

			Assembly a0 = (Assembly) assems[0];
			Assembly a1 = (Assembly) assems[1];
			selfBSs0 = a0.getBSCounts();
			selfBSs1 = a1.getBSCounts();

			updateCounts(otherFreeBSs, selfBSs0, false);
			updateCounts(otherFreeBSs, selfBSs1, false);

			if (a0.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(otherMonomerBSs, selfBSs0, false);
			}
			if (a1.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(otherMonomerBSs, selfBSs1, false);
			}

			a0.setValidTime(curtime);
			a1.setValidTime(curtime);
			Event f0 = null;
			Event f1 = null;

			f0 = getFormBondEvent(a0, otherFreeBSs, otherMonomerBSs, a1);
			f1 = getFormBondEvent(a1, otherFreeBSs, otherMonomerBSs, a0);

			Event c = sampleTwoAssems(a0, a1);
			myMinEvent0 = screenEvents(f0, c);
			myMinEvent0 = screenEvents(myMinEvent0, getConfChangeEvent(a0));
			Event tbe0 = sampleBreakBondEvent(a0);
			myMinEvent0 = screenEvents(myMinEvent0, tbe0);

			if (a0.numSubunits() == Test.sizeOfSubunit) {
				if (tbe0 != null)
					System.out.println(a0.numSubunits() + " " + currentStep
							+ " errorin2simulation");
			}

			myMinEvent1 = screenEvents(f1, c);
			myMinEvent1 = screenEvents(myMinEvent1, getConfChangeEvent(a1));
			Event tbe1 = sampleBreakBondEvent(a1);
			myMinEvent1 = screenEvents(myMinEvent1, tbe1);

			if (a1.numSubunits() == Test.sizeOfSubunit) {
				if (tbe1 != null)
					System.out.println("errorin3simulation");
			}

		}

		if ((myMinEvent0 != null) && (myMinEvent1 == null)) {
			sendEvent(myMinEvent0);
		} else if ((myMinEvent0 == null) && (myMinEvent1 != null)) {
			sendEvent(myMinEvent1);
		} else if ((myMinEvent0 != null) && (myMinEvent1 != null)) {
			if (myMinEvent0 != myMinEvent1) {
				sendEvent(myMinEvent0);
				sendEvent(myMinEvent1);
			} else {
				sendEvent(myMinEvent0);
			}
		} else {
		}
	}



	/**
	 * Using the two assemblies given to picks a new FormBondEvent. This
	 * will be the minimum of the possible events for these Assemblies.
	 * 
	 * @param assem - The fist Assembly to sample
	 * @param otherFreeBSs - A count of the free BindingSites 
	 * 		for the first Assembly
	 * @param otherMonomerBSs - A count of free BindingSites 
	 * 		for the second Assembly
	 * @param assem2 - The second Assembly to sample
	 * @return FormBondEvent - the new minimum bonding event.
	 */
	private Event getFormBondEvent(Assembly assem,
			HashMap otherFreeBSs, HashMap otherMonomerBSs, Assembly assem2) {

		double duration = Double.MAX_VALUE;
		double tem;
		String bstSelfChoice = null;
		String bstOtherChoice = null;
		int countSelf = -1;
		int countOther = -1;

		HashMap<String, Integer> abs = assem.getBSCounts();
		HashMap<String, Integer> other = null;

		boolean monomerOnly = false;

		if ((assem.numSubunits() == Test.sizeOfSubunit) || (!Test.bindMonomerOnly)) 
			other = otherFreeBSs;
		else {
			//current Assembly has more than one monomers and
			// Test.bindMonomerOnly=true
			monomerOnly = true;
			other = otherMonomerBSs;
		}

		Iterator<String> it = abs.keySet().iterator();
		while (it.hasNext()) {

			String bst = it.next();
			BindingSiteType bsType = bstMap.get(bst);
			Iterator<String> partner = bsType.getPartners().iterator();
			while(partner.hasNext())
			{

				String bstP = partner.next();

				int a = abs.get(bst).intValue();
				int b = other.get(bstP).intValue();
				int product = a * b;

				if (product != 0) {

					double bindTime = bsType.getBindTime(bstP);
					tem = getExp(bindTime/product);

					if (tem < duration) {
						duration = tem;
						countSelf = a;
						countOther = b;
						bstSelfChoice = bst;
						bstOtherChoice = bstP;
					}
				}
			}
		}

		Event minFormBondEvent = null;

		if (duration < Double.MAX_VALUE) {
			minFormBondEvent = sampleFormBondEvent(duration, assem, assem2,
					bstSelfChoice, bstOtherChoice,
					(int) (countSelf * getRand()),
					(int) (countOther * getRand()), monomerOnly);
		}
		return minFormBondEvent;
	}

	/**
	 * A Simple method to screen to <code>Event</code>s. Will return the
	 * second if the first is null, first if second is null, or the minimum
	 * of the two if both are not null
	 * 
	 * @param e0 - The first <code>Event</code>
	 * @param e1 - The second <code>Event</code>
	 * @return <code>Event</code> - Whichever <code>Event</code> 
	 * 	is not null and occurs first. 
	 */
	private Event screenEvents(Event e0, Event e1) {

		if (e0 == null) 
			return e1;
		else if (e1 == null) 
			return e0;
		else if (e0.compareTo(e1) == 1) 
			return e0;
		else
			return e1;
	}



	/**
	 * Samples all possible <code>BreakBondEvents</code> for a given
	 * <code>Assembly</code>
	 * 
	 * @param assem - The <code>Assembly</code> to sample
	 * @return <code>BreakBondEvents</code> if an event is found, <code>null</code>
	 * 	otherwise
	 */
	private Event sampleBreakBondEvent(Assembly assem) {

		Vector<Subunit> subs = assem.getSubunits();

		if (subs.size() == Test.sizeOfSubunit) 
			return null;

		double duration = 0;
		double minTime = Double.MAX_VALUE;

		BindingSite temBs, temPartner;
		BindingSite minBs = null;
		BindingSite minPartner = null;
		CheckPairSeen cps = new CheckPairSeen();

		//go thru boundsubs, looking at all possible bond breaking events
		int size = subs.size();
		for (int j = 0; j < size; ++j) {

			Subunit sub = subs.get(j);
			Vector<BindingSite> bindsites = sub.getBindingSites();

			if ((sub.getBoundSubunits().size() > 1) && (Test.breakOnlyEnds)) {
				continue;
			}
			int bsSize = bindsites.size();
			for (int k = 0; k < bsSize; ++k) {

				temBs = bindsites.get(k);

				if (temBs.isBound()) {

					temPartner = temBs.getPartner();

					if (cps.marked(temBs.getID(), temPartner.getID()))
						continue;

					//else mark it seen now
					cps.addPair(temBs.getID(), temPartner.getID());

					if (Test.noLoopOnly) {
						boolean breakLoop = assem.splitAssemblyInLoop(temBs,
								temPartner);
						if (breakLoop)
							continue;
					}

					double brkTime =bstMap.get(temBs.getBSTName()).getBreakTime(temPartner.getBSTName());

					duration = getExp(brkTime);
					//if this is the fastest Event so far, store it
					if (minTime >= duration) {
						minTime = duration;
						minBs = temBs;
						minPartner = temPartner;
					}
				}
			}
		}


		if (minTime == Double.MAX_VALUE)
			return null;

		Assembly[] assems = { minBs.getAssembly(), null };

		return new Event(curtime, minTime + curtime, assems, minBs,
				minPartner,EventType.brkBndEvt);
	}



	private Event sampleFormBondEvent(double duration, Assembly assem,
			Assembly assem2, String bstSelf, String bstOther, int indexSelf,
			int indexOther, boolean monomerOnly) {

		/**
		 * 
		 * Walks through all Assembly(s) in this Simulation, finds and saves all
		 * the other free BindingSite(s) by checking if each one is bound and is
		 * in Assembly(s) specified by assems
		 */
		BindingSite minBs = null;
		BindingSite minPartner = null;
		Vector<Subunit> subs = assem.getSubunits();
		boolean indexSelfFound = false;
		int temIndex = 0;
		int aSize = subs.size();
		for (int i = 0; ((i < aSize) && (!indexSelfFound)); ++i) {

			Subunit sub = subs.get(i);
			Vector<BindingSite> bss = sub.getBindingSites();
			int bSize = bss.size();
			for (int j = 0; j < bSize && !indexSelfFound; ++j) {

				BindingSite bs = bss.get(j);

				if ((!bs.isBound()) && (bs.getBSTName().equals(bstSelf))) {

					if (temIndex == indexSelf) {
						indexSelfFound = true;
						minBs = bs;
					} else {
						temIndex++;
					}
				}
			}
		}

		//find partner BindingSite of type bstOther from all the other
		// Assembly(s) or all the other monomers if argument monomerOnly is true
		boolean indexOtherFound = false;
		temIndex = 0;

		Iterator<Assembly> it = assembliesHashMap.values().iterator();
		Assembly temAssembly;

		while (it.hasNext() && !indexOtherFound) {

			temAssembly = it.next();

			if (temAssembly == assem || temAssembly == assem2)
				continue;
			if (temAssembly.numSubunits() > Test.sizeOfSubunit && monomerOnly) 
				continue;

			Vector<Subunit> temSubs = temAssembly.getSubunits();
			int cSize = temSubs.size();
			for (int j = 0; ((j < cSize) && (!indexOtherFound)); ++j) {

				Subunit temSub = (Subunit) temSubs.get(j);
				Vector<BindingSite> temBss = temSub.getBindingSites();
				int dSize = temBss.size();
				for (int k = 0; ((k < dSize) && (!indexOtherFound)); ++k) 
				{
					BindingSite temBs = temBss.get(k);

					if ((!temBs.isBound()) && (temBs.getBSTName().equals(bstOther))) {

						if (temIndex == indexOther) {
							indexOtherFound = true;
							minPartner = temBs;							
						} else 
							temIndex++;
					}
				}
			}
		}

		Assembly[] assems = { minBs.getAssembly(), 
				minPartner.getAssembly() };

		return new Event(curtime, curtime + duration, assems, minBs,
				minPartner,EventType.formBndEvt);
	}

	public final static boolean isCompatible(String bst1, String bst2)
	{
		return bstMap.get(bst1).isCompatible(bst2);
	}


	private Event sampleTwoAssems(Assembly A1, Assembly A2) {
		//Assume A1 and A2 should be different Assembly(s)
		if (A1.equals(A2)) {
			return null;
		}

		if ((A1.numSubunits() > Test.sizeOfSubunit)
				&& (A2.numSubunits() > Test.sizeOfSubunit)
				&& (Test.bindMonomerOnly))
			return null;

		Vector<BindingSite> bss1 = A1.getFreeBindingSites();
		Vector<BindingSite> bss2 = A2.getFreeBindingSites();

		int primary = -1;
		int partner = -1;

		double minTime = Double.MAX_VALUE;
		double duration = 0;
		int aSize = bss1.size();
		for (int i = 0; i < aSize; ++i) {
			BindingSite temBs = bss1.get(i);
			int bSize = bss2.size();
			for (int j = 0; j < bSize; ++j) {
				BindingSite temPartner = bss2.get(j);
				//check for compatible bindingsiteTypes

				if (!isCompatible(temBs.getBSTName(), temPartner.getBSTName())) {
					continue;
				}


				double bTime = bstMap.get(temBs.getBSTName()).getBindTime(temPartner.getBSTName());


				duration = getExp(bTime);

				//if this is the fastest Event so far, store it
				if (minTime >= duration) {
					minTime = duration;
					primary = i;
					partner = j;
				}
			}
		}
		//To be removed for new initialize()
		Assembly[] assemsInvolved = { A1, A2 };

		if (minTime != Double.MAX_VALUE) {
			return new Event(curtime, curtime + minTime,
					assemsInvolved, (BindingSite) bss1.get(primary),
					(BindingSite) bss2.get(partner),EventType.formBndEvt);
		} else {
			return null;
		}

	}



	/**
	 * This will pick a <code>ConfChangeEvent</code> for the given
	 * <code>Assembly</code>. Assembly should be a monomer, as defined 
	 * by Test.sizeofSubunit.
	 * 
	 * @param assem_i - The <code>Assembly</code> to pick a <code>Event</code>
	 * 	for.
	 * @return <code>ConfChangeEvent</code> if possible, null if none
	 * 	found.
	 */
	private Event getConfChangeEvent(Assembly assem_i) {

		if (!Test.csAllowed || assem_i.numSubunits() != Test.sizeOfSubunit)
			return null;

		Vector<Subunit> subs = assem_i.getSubunits();
		Subunit sub = subs.get(0);

		Vector<Domain> domains = sub.getDomains();
		double minTime = Double.MAX_VALUE;
		int dID = -1;
		Conformation conf = null;

		/**
		 * Go through all domains in the Subunit and look at all possible
		 * conformations
		 */
		int aSize = domains.size();
		for (int i = 0; i < aSize; ++i) {

			Domain curDomain = domains.get(i);
			int id = curDomain.getDomainId();
			Conformation curConf = curDomain.getCurConf();
			HashMap<String,Conformation> allConfs = curDomain.getConfs();
			HashMap<String, Double> switchMap = confTimes.get(curConf.getName());

			/**
			 * For each possible conformation of each possible domain use the
			 * switch map to find the time for the event, storing the smallest
			 */
			Iterator<Conformation> cnfItr = allConfs.values().iterator();
			while(cnfItr.hasNext())
			{
				Conformation testConf = cnfItr.next();

				if (switchMap.containsKey(testConf.getName())) {

					double changeTime = switchMap.get(testConf.getName()).doubleValue();	
					double duration = getExp(changeTime);

					if (duration < minTime) {
						minTime = duration;
						dID = id;
						conf = testConf;
					}
				}
			}
		}

		if (minTime == Double.MAX_VALUE)
			return null;

		return new Event(curtime, curtime + minTime,
				new Assembly[] { assem_i, null}, sub,
				dID, conf);
	}




	private void updateCounts(HashMap<String, Integer> old, HashMap<String, Integer> tem, boolean add) {

		Iterator<String> it = tem.keySet().iterator();

		while (it.hasNext()) {

			String name = it.next();

			int temN = tem.get(name).intValue();
			int oldN = old.get(name).intValue();

			if (add) 
				old.put(name, new Integer(temN + oldN));
			else 
				old.put(name, new Integer(oldN - temN));
		}
	}



	private void simpleUpdateCountsAll(BindingSite b, boolean add) {

		String name = b.getBSTName();
		int n = freeBSs.get(name).intValue();

		if (add)
			freeBSs.put(name, new Integer(n + 1));
		else
			freeBSs.put(name, new Integer(n - 1));
	}



	/**
	 * Returns an random exponentially distributed variable with the given
	 * average value.
	 * 
	 * @param avg - The average to base this around
	 * @return double - a random exponential
	 */
	private double getExp(double avg) {

		double t = -avg * Math.log(getRand());
		return t;
	}



	/**
	 * Returns a random double within [0 , 1)
	 * 
	 * @return double
	 */
	private double getRand() 
	{ return Test.rand.nextDouble(); }





	/**
	 * Adds an <code>Event</code> to the queue
	 * 
	 * @param e - The <code>Event</code> to add
	 */
	private final void sendEvent(Event e) 
	{ pq.add(e); }


	/**
	 * This method writes the size distribution. It will write to 
	 * file if <code>Test.printToScreen</code> = false, otherwise it will print
	 * to screen. The output is formatted as follows:
	 * 		Simulation_Time Count_of_Assemblies_of_Size_1 ... Count_of_Assemblies_of_Size_N
	 * where N is defined by <code>Test.maxOutputSize</code>
	 */
	public void printDistributionToScreen() {

		StringBuffer output = new StringBuffer();

		output.append(curtime + " ");

		for (int i = 0; i < currentSizeDistribution.length && 
		i < Test.maxOutputSize; ++i)
		{
			output.append(currentSizeDistribution[i] + " ");
			//int count = currentSizeDistribution[i];
			//if(count!=0)
			//	output.append(" ["+(i+1)+"]"+count + " ");
		}

		output.append("\n");
		System.out.print(output);
	}

	private void printToFile()
	{
		if( (got73 && got80) || fileCounter >3) return;
		
		Iterator<Assembly> aItr = assembliesHashMap.values().iterator();
		
		while(aItr.hasNext())
		{
			Assembly a  = aItr.next();
			if(a.numSubunits()>= 73)
			{
				if(a.numSubunits() >=73 && a.numSubunits()  < 80){
					if(got73) return;
					got73 = true;					
				}
				
				if(a.numSubunits() >=80 ){
					if(got80) return;
					got80 = true;
				}
				FileWriter fstream = null;
				try {
					File f = new File(fileCounter.toString().concat(".txt"));
					fstream = new FileWriter(f);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				BufferedWriter out = new BufferedWriter(fstream);
				StringBuilder sb = new StringBuilder();
				Iterator<Subunit> subs = a.getSubunits().iterator();
				while(subs.hasNext())
				{
					Subunit s = subs.next();
					Iterator<BindingSite> bss = s.getBindingSites().iterator();

					while(bss.hasNext())
					{
						BindingSite b1 = bss.next();
						sb.append(s.getID());sb.append(",");
						sb.append(s.getPositionReal().toString().replace("(", "").replace(")",""));sb.append(",");
						sb.append(b1.getID());sb.append(",");
						sb.append(b1.getPosReal().toString().replace("(", "").replace(")",""));sb.append(",");
						if(!typemap.containsKey(b1.getBSTName()))
						{
							System.out.println("Error in: "+b1.getBSTName());
							System.exit(0);
						}
						sb.append(typemap.get(b1.getBSTName()));
						sb.append(",");
						if(b1.isBound())
						{
							BindingSite b2 = b1.getPartner();
							Subunit s2 = b2.getSubunit();
							sb.append(s2.getID());sb.append(",");
							sb.append(s2.getPositionReal().toString().replace("(", "").replace(")",""));sb.append(",");
							sb.append(b2.getID());sb.append(",");
							sb.append(b2.getPosReal().toString().replace("(", "").replace(")",""));sb.append(",");
							if(!typemap.containsKey(b2.getBSTName()))
							{
								System.out.println("Error in: "+b2.getBSTName());
								System.exit(0);
							}
							sb.append(typemap.get(b2.getBSTName()));
														
						}
						else{
							sb.append("inf,inf,inf,inf,inf,inf,inf,inf,inf");
						}
						sb.append("\n");
					}
				}
				try {
					
					out.write(sb.toString());
					out.close();
					++fileCounter;
					/*if(fileCounter>3)
						System.exit(-1);*/
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * A class used when selecting <code>BreakBondEvents</code> to store what
	 * <code>BindingSites</code> have already been seen. This is done by storing
	 * their ID.
	 * 
	 * @author Tiequan Zhang
	 * @author Blake Sweeney
	 * @version 1.4
	 */
	private class CheckPairSeen {

		/** vector of double[2], which are the seen pairs */
		private Vector<double[]> seen; 

		/**
		 * A Constructor
		 */
		public CheckPairSeen() {
			seen = new Vector<double[]>();
		}



		/**
		 * Adds a seen pair
		 * 
		 * @param i - first ID
		 * @param j - second ID
		 */
		public void addPair(int i, int j) { 
			seen.add(new double[] { i, j });
		}

		/**
		 * Checks if a pair has been seen. Both must have been seen together
		 * for it to be considered seen.
		 * 
		 * @param i - First ID to look for
		 * @param j - Second ID to look for
		 * @return <code>boolean</code> - true if both ID's have been seen together
		 */
		public boolean marked(int i, int j) {

			double[] tmp;
			int s = seen.size();
			for (int k = 0; k < s; ++k) {
				tmp = (double[]) (seen.get(k));
				if ((tmp[0] == i && tmp[1] == j)
						|| (tmp[0] == j && tmp[1] == i)) {
					return true;

				}
			}
			return false;
		}


	}
}
