package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

public class StsSensorVirtualVolumeClass extends StsVirtualVolumeClass implements StsClassTimeDisplayable, StsSerializable, StsClassSurfaceDisplayable, StsClassCursor3dTextureDisplayable //, StsClassDisplayable
{
	public StsSensorVirtualVolumeClass()
	{
		userName = "Virtual Volume from Sensors";
	}

	public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
				{
						new StsBooleanFieldBean(this, "isVisibleOnCursor", "On 3D Cursors"),
						new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
						new StsBooleanFieldBean(this, "contourColors", "Contoured Seismic Colors")
				};
	}

	public void setIsVisibleOnCursor(boolean visible)
	{
		if (isVisibleOnCursor == visible)
			return;
		isVisibleOnCursor = visible;
		currentModel.win3dDisplayAll();
	}

	public boolean getIsVisibleOnCursor()
	{
		return isVisibleOnCursor;
	}

	public void setContourColors(boolean contour)
	{
		if (this.contourColors == contour)
			return;
		this.contourColors = contour;
		currentModel.win3dDisplayAll();
	}

	public boolean getContourColors()
	{
		return contourColors;
	}

	public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
	{
		StsObject[] objects = getObjectList();
		for (int i = 0; i < objects.length; i++)
		{
			((StsSensorVirtualVolume) objects[i]).recomputeSensorRanges();
		}
	}
}
