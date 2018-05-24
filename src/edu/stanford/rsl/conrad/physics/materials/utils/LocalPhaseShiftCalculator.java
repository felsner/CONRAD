package edu.stanford.rsl.conrad.physics.materials.utils;


//similar to package edu.stanford.rsl.conrad.physics.materials.utils.LocalMassAttenuationCalculator

//created by Lina Felsner 22.02.2018

public class LocalPhaseShiftCalculator {
	

	/**
	 * Reingehackt. 
	 * 
	 * @param comp is the {@link WeightedAtomicComposition} of the material
	 * @param energy is the energy of interest
	 * @param attType is the {@link AttenuationType} of interest
	 * @return 
	 */
	public static double getDeltaValues(String name, WeightedAtomicComposition comp, double energy, AttenuationType attType, double density) {
		
		// abtesten welcher compound
		double delta = 0;
		
		
		if(energy == 0.082) {
			// Werte ensprechen delta_theo mit 82keV
			// Quelle: Willner Paper 'Quantitative X-ray phase-contrast computed tomography at 82keV'
			switch (name.toLowerCase()) {
			case "water":
				delta = 3.96E-08;
				// TODO
				//System.out.println("fix me! instead of the delta value for water PMMA is used.");
				break;
		    case "pvc":
		    	delta = 4.43E-08;
		    	break;
		    case "ptfe":
		    	delta = 6.52E-08;
		    	break;
		    default: 
		    	//System.out.println("warning! phase values are not known. Use Attenuation values...");
		    	delta = density * LocalMassAttenuationCalculator.getMassAttenuationData(comp, energy, attType);
		    	break;
			}
		}else if(energy == 0.051){
			// Werte ensprechen delta_theo mit 51keV
			// Quelle: Sarapata Paper 'Quantitaive imaging using high-energy X-ray phase-contrast CT with a 70 kVp polychromatic x-ray spectrum'
			switch (name.toLowerCase()) {

		    case "pvc":
		    	delta = 1.08E-07;
		    	break;
		    case "ptfe":
		    	delta = 1.60E-07;
		    	break;
		    case "silicondioxid":
				delta = 1.66E-07;
				break;	
		    default: 
		    	//System.out.println("warning! phase values are not known. Use Attenuation values...");
		    	delta = density * LocalMassAttenuationCalculator.getMassAttenuationData(comp, energy, attType);
		    	break;
			}
			
		} else {
			delta = density * LocalMassAttenuationCalculator.getMassAttenuationData(comp, energy, attType);
		}
		
		return delta;
	}



}
