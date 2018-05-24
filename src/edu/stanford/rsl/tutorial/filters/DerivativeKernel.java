package edu.stanford.rsl.tutorial.filters;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * First part of the decomposition of the ramp filter into a first derivative and a Hilbert filter.
 * See L. Zeng. "Medical Image Reconstruction: A Conceptual tutorial". 2009, page 28
 * @author akmaier
 *
 */

public class DerivativeKernel implements GridKernel {
	
	private Grid1D dKernel;
	public DerivativeKernel() {
		dKernel = new Grid1D(new float[3]);
		
		dKernel.setAtIndex(0, 0.f);// fix by Lina Felsner. Convention: add sample at the beginning
		dKernel.setAtIndex(1, -1.f);
		dKernel.setAtIndex(2, 1.f); 
		// convolveFloat expects an odd kernel size
		
	}

	public void applyToGrid(Grid1D input) {
		float[] inputFloat = input.getBuffer();
		ImageProcessor ip = new FloatProcessor(inputFloat.length, 1, inputFloat);
		Convolver c = new Convolver();
		c.convolveFloat((FloatProcessor) ip, dKernel.getBuffer(), 3, 1);
		NumericPointwiseOperators.copy(input, new Grid1D(inputFloat));
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Grid1D data = new Grid1D(new float[100]);
		
		for (int i=0;i<data.getSize()[0];++i) data.setAtIndex(i, (float) i);
		DerivativeKernel dKern = new DerivativeKernel();
		dKern.applyToGrid(data);
		VisualizationUtil.createPlot(data.getBuffer()).show();
	}

}
/*
 * Copyright (C) 2010-2014 Andreas Maier
 * CONRAD is developed as an Open Source project under the GNU General Public License (GPL).
*/