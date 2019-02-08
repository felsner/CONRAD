/*
 * Copyright (C) 2014  Shiyang Hu
 * CONRAD is developed as an Open Source project under the GNU General Public License (GPL).
 */
package edu.stanford.rsl.conrad.data.numeric;


/**
 * Grid class to model multi channel data. It can be used to model vector valued voxels, but will still fit into
 * the image streaming concept.
 * 
 * @author akmaier
 *
 */
public class MultiChannelGrid3D extends Grid3D {

	String [] channelNames;
	Grid4D multichannelData;
	
	public MultiChannelGrid3D(int width, int height, int depth, int channels) {
		// This will cause some parent methods to fail ( e.g. getSize() );
		// TODO Override necessary parent methods.
		super(0,0,0);
		multichannelData = new Grid4D(width, height, depth, channels);
		buffer = multichannelData.getSubGrid(0).getBuffer();
		//this(width, height, depth, true);

		// Lina
		this.size = new int[] { width, height, depth};
		this.spacing = new double[] {1,1,1};
		this.origin = new double[] {0,0,0};
	}
	
	/**
	 * Returns the Grid3D of the respective Channel
	 * @param c the channel number
	 * @return the Grid3D
	 */
	public Grid3D getChannel(int c){
		Grid3D intermediate = multichannelData.getSubGrid(c);
		intermediate.setSpacing(this.getSpacing());
		intermediate.setOrigin(this.getOrigin());
		return intermediate;
	}
	
	/**
	 * Sets the channel in the image to the given grid 2D
	 * @param c the channel number
	 * @param channel the grid 2D
	 */
	public void setChannel(int c, Grid3D channel){
		multichannelData.setSubGrid(c, channel);
	}
	
	/**
	 * Reports the number of channels
	 * @return the number of channels
	 */
	public int getNumberOfChannels(){
		return multichannelData.getSize()[3];
	}
	
	
	
	/**
	 * Set a pixel value at position (x,y)
	 * @param x The value's x position
	 * @param y The value's y position
	 * @param c The value's channel
	 * @param value The value to set
	 */
	public void putPixelValue(int x, int y, int z, int c, float value) {
		multichannelData.setAtIndex(x, y, z, c, value);
	}
	
	
	/**
	 * Set a pixel value at position (x,y)
	 * @param x The value's x position
	 * @param y The value's y position
	 * @param c The value's channel
	 * @param value The value to set
	 */
	public void putPixelValue(int x, int y, int z, int c, double value) {
		multichannelData.setAtIndex(x, y, z, c, (float)value);
	}
	
	
	/**
	 * Get the pixel value at position (x,y)
	 * @param x The value's x position
	 * @param y The value's y position
	 * @param c The value's channel
	 * @return the value of the pixel
	 */
	public float getPixelValue(int x, int y, int z, int c) {
		return multichannelData.getAtIndex(x, y, z, c);
	}

	@Override
	public NumericGrid clone() {
		return (new Grid3D(this));
	}

	/**
	 * @return the channelNames
	 */
	public String[] getChannelNames() {
		return channelNames;
	}

	/**
	 * @param channelNames the channelNames to set
	 */
	public void setChannelNames(String[] channelNames) {
		this.channelNames = channelNames;
	}

	@Override
	/**
	 * @return The array's size (excluding borders).
	 */
	public int[] getSize() {
		return new int[]{this.multichannelData.getSize()[0],this.multichannelData.getSize()[1],this.multichannelData.getSize()[2]};
	}
	
	
	@Override
	public void setSpacing(double... spacing){
		super.setSpacing(spacing);
		multichannelData.setSpacing(spacing[0], spacing[1], spacing[2], 0.0);
	}
	
	@Override
	public void setOrigin(double... origin){
		super.setOrigin(origin);
		multichannelData.setOrigin(origin[0], origin[1], origin[2], 0.0);
	}
}
