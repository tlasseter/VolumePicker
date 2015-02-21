package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SurfaceCurvature.*;
import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.nio.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 14, 2009
 * Time: 10:35:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPatchVolume extends StsSeismicBoundingBox implements StsTreeObjectI, StsSerializable
{
	/** patches created are sorted by row, col, and z in this array */
	public StsPatchGrid[] rowSortedPatchGrids;
	protected StsColorscale colorscale;
	protected StsHistogram histogram;
	protected StsSeismicVolume seismicVolume;
	protected String seismicName;
	protected boolean filter = false;
	protected int boxFilterWidth = 1;
	transient public int nInterpolatedSlices;
	transient public float interpolatedZInc;
	transient public int interpolatedSliceMin;
	transient public int interpolatedSliceMax;
	transient StsPatchVolumeClass patchVolumeClass;
	transient public StsCroppedBoundingBox croppedBoundingBox;

	/**
	 * gridList contains patches which have been completed. At the end of a row, prevRowGrids contains disconnected grids which
	 * are added to gridList patches.  At the end of all rows, remaining grids are in rowGrids which are then added to gridList.
	 */
	transient ArrayList<StsPatchGrid> gridList;

	/**
	 * rowGrids contains new patches and existing patches from previous row connected to the latest row;
	 * at the completion of the row, these become the previousRowGrids. At start of row, its initialized to empty.  If a new grid is created,
	 * it is added to rowGrids.  If an existing grid is connected to a window in the row, it is added to rowGrid and removed from prevRowGrid.
	 * At the end of the row, grids still connected are in rowGrid, and disconnected ones are in prevRowGrids. These disconnected grids are
	 * added to gridList.  prevRowGrids is then set to rowGrids and rowGrids initialized for the next row.
	 */
	transient HashMap<Integer, StsPatchGrid> rowGrids = null;
	transient Iterator<StsPatchGrid> rowGridsIterator;

	/** prevRowGrids are the active patches in the previousRow; when making a connection to a window in the previous row, we look here for a patch. */
	transient HashMap<Integer, StsPatchGrid> prevRowGrids = null;
	transient Iterator<StsPatchGrid> prevRowGridsIterator;

	/** a temporary array built at initialization with final patches sorted by col, row, and z */
	transient public StsPatchGrid[] colSortedPatchGrids;
	/** total number of points on all patches on this volume */
	transient int nPointsTotal;
	/** window size in wavelengths. If zero, window size is 1/2 wavelength. */
	transient int windowSize;
	/** half window size in wavelengths. */
	transient float halfWindowSize;
	/**
	 * Minimum data amplitude that will be used for a minimum or maximum event.
	 * Max data amplitude is Min(-dataMin, dataMax). So to be used, the abs(amplitude) >= minDataAmplitude
	 */
	transient double minDataAmplitude;
	/** pick window on adjoining trace cannot be more than this many wavelengths away from picked window */
	transient float pickDifWavelengths;
	/** Type(s) of pick events to be correlated: min, max, min&max, all */
	transient byte pickType;
	/** number of points to find in half of window */
	transient int nHalfSamples;
	/** Window center is either max or min; window ends are either the same or a zero crossing;  This flag indicates window ends with zero-crossing */
	transient boolean windowEndIsZeroCrossing;
	/** multiplier of half window size yielding the max allowed pick difference */
	transient double halfWindowPickDifFactor;
	/** number of interpolation intervals between samples */
	transient int nInterpolationIntervals;
	/** indicates iterative stretchCorrel is to be applied from max down to min by inc */
	transient boolean isIterative;
	/** max value for stretchCorrelation */
	transient float autoCorMax;
	/** min value for stretchCorrelation */
	transient float autoCorMin;
	/** increment of stretch correlation in interative pick */
	transient float autoCorInc;
	/** manual picking minimum acceptable cross-correlation */
	transient float manualCorMin;
	/** sequence of stretchCorrelations: 1 if not iterative, max to min by inc if iterative */
	transient float[] stretchCorrelations;
	/** number of stretchCorrelations in sequence */
	transient int nIterations;
	/** correlate using falseTypes (e.g., a false Max matches a Max, etc) */
	transient boolean useFalseTypes = false;
	/** double check connection by matching it back from selected otherWindow; accept match if backMatchWindow is null or the same */
	transient boolean checkBackMatch = true;
	/** For debugging: only run the cycle skip check if true */
	transient boolean backMatch = true;
	/** row currently being computed: used for debugPatchGrid print out only */
	transient int row, col, volRow, volCol;
	transient int nPatchPointsMin;
	transient public byte curveType = CURVPos;
	transient public int filterSize = 0;
	transient protected boolean displaySurfs = false;
	transient protected boolean displayVoxels = false;
	transient int nSmallGridsRemoved;
	transient int nParentGrids;
	transient StsPoint currentCursorPoint;
	transient StsPatchGrid cursorPointPatch;
	transient private StsPatchGrid[] selectedPatchGrids;

	transient float[] histogramValues;
	transient int nHistogramValues = 10000;

	/** if the row or col correl is < this value, that link and the associated edge is not drawn */
	static public float minLinkCorrel = 0.0f;

	static protected StsObjectPanel objectPanel = null;

	public static final byte PICK_MAX = 0;
	public static final byte PICK_MIN = 1;
	public static final byte PICK_MIN_MAX = 2;
	public static final byte PICK_ALL = 3;

	public static final int largeInt = 99999999;

	public static final String[] pickTypeNames = new String[]{"All", "Min+Max", "Maximum", "Minimum"}; //, "Zero-crossing+", "Zero-crossing-", "All"};
	public static final String[] stopCriteriaNames = new String[]{"Stop", "Replace", "Stop if same Z"};
	//static StsComboBoxFieldBean displayAttributeBean = new StsComboBoxFieldBean(StsPatchVolume.class, "displayAttribute", "Attribute");
	static public final float badCurvature = StsPatchGrid.badCurvature;

	public String getSeismicName()
	{
		return seismicName;
	}

	public void setSeismicName(String seismicName)
	{
		this.seismicName = seismicName;
	}

	static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsPatchVolume.class, "colorscale");


	static public final StsFieldBean[] displayFields =
			{
					new StsBooleanFieldBean(StsPatchVolume.class, "isVisible", "Display on Cursors"),
					new StsBooleanFieldBean(StsPatchVolume.class, "displaySurfs", "Display as Surfaces"),
					new StsBooleanFieldBean(StsPatchVolume.class, "displayVoxels", "Display as Voxel cloud"),
			};

	static public final StsFieldBean[] propertyFields = new StsFieldBean[]
			{
					new StsStringFieldBean(StsPatchVolume.class, "name", true, "Name"),
					colorscaleBean
			};

	final public float getZ(int slice)
	{
		return zMin + this.interpolatedZInc * slice;
	}

	public float getDataMin()
	{
		return dataMin;
	}

	public float getDataMax()
	{
		return dataMax;
	}

	static public final int defaultTestWindow = 21;
	static public final float defaultTestWavelengths = 1.0f;

	//Curvature Attribute Types
	static public final byte CURVDip = StsSurfaceCurvatureAttribute.CURVDip;
	static public final byte CURVStrike = StsSurfaceCurvatureAttribute.CURVStrike;
	static public final byte CURVMean = StsSurfaceCurvatureAttribute.CURVMean;
	static public final byte CURVGauss = StsSurfaceCurvatureAttribute.CURVGauss;
	static public final byte CURVPos = StsSurfaceCurvatureAttribute.CURVPos;
	static public final byte CURVNeg = StsSurfaceCurvatureAttribute.CURVNeg;
	static public final byte CURVMin = StsSurfaceCurvatureAttribute.CURVMin;
	static public final byte CURVMax = StsSurfaceCurvatureAttribute.CURVMax;


	static final float largeFloat = StsParameters.largeFloat;

	/** debugPatchGrid prints showing row operations */
	static final boolean debug = true;
	/** turn on timer for curvature operation */
	static final boolean runTimer = false;
	/** millisecond timer */
	static StsTimer timer;
	/** debugPatchGrid for tracePoints linkList operations */
	static final boolean debugTracePointsLink = false;
	/** debugPatchGrid for CorrelationWindows */
	static final boolean debugCorrelationWindows = false;
	/** search for multiple windows */
	static final boolean searchForMultipleWindowMatches = true;
	/** print patch operations and draw only this patch */
	static boolean drawPatchBold = debug && StsPatchGrid.debugPatchGrid; // StsPatchGrid.debugPatchGrid;
	/** debug: point clone operations */
	static final boolean debugCloneOK = true;
	/** debug: connect closest points only */
	static final boolean debugConnectCloseOnly = true;
	static public String iterLabel = "";

	private static final long serialVersionUID = 1L;

	public StsPatchVolume()
	{
	}

	public StsPatchVolume(StsSeismicVolume seismicVolume)
	{
		super(false);
		StsToolkit.copySubToSuperclass(seismicVolume, this, StsRotatedGridBoundingBox.class, StsBoundingBox.class, true);
		this.seismicVolume = seismicVolume;
		seismicName = seismicVolume.getName();
		zDomain = seismicVolume.zDomain;
		stsDirectory = seismicVolume.stsDirectory;
		initialize(currentModel);
	}

	/**
	 * This is central method for constructing the volume of patchGrids.
	 * For each row, we examine correlation with traces in same row & prev col and same col & prev row.
	 * @param pickPanel graphics panel with progress bar updated as each row completed
	 */
	public void constructPatchVolume(StsPatchPickPanel pickPanel)
	{
		StsProgressBar progressPanel = pickPanel.progressBar;
		windowSize = pickPanel.corWavelength;
		pickDifWavelengths = pickPanel.maxPickDif;
		float minAmpFraction = pickPanel.minAmpFraction;
		float maxStretch = pickPanel.maxStretch;
		pickType = pickPanel.pickType;
		nPatchPointsMin = pickPanel.minPatchSize;
		useFalseTypes = pickPanel.useFalseTypes;
		checkBackMatch = pickPanel.checkBackMatch;
		backMatch = pickPanel.checkBackMatch;
		isIterative = pickPanel.isIterative;
		autoCorMax = pickPanel.autoCorMax;
		autoCorMin = pickPanel.autoCorMin;
		autoCorInc = pickPanel.autoCorInc;
		manualCorMin = pickPanel.manualCorMin;

		StsPatchGrid.initializeDebug(pickPanel);
		initializeLists();

		if (!isIterative)
		{
			nIterations = 1;
			stretchCorrelations = new float[]{manualCorMin};
		}
		else
		{
			nIterations = StsMath.ceiling(1 + (autoCorMax - autoCorMin) / autoCorInc);
			stretchCorrelations = new float[nIterations];
			float stretchCorrelation = autoCorMax;
			for (int n = 0; n < nIterations; n++, stretchCorrelation -= autoCorInc)
				stretchCorrelations[n] = stretchCorrelation;
		}

		rowSortedPatchGrids = new StsPatchGrid[0];
		colSortedPatchGrids = null;
		initialize();
		StsPatchGrid.staticInitialize();

		if (progressPanel != null)
			progressPanel.initialize(croppedBoundingBox.nRows);

		TracePoints[] rowTracePoints = null; //new TracePoints[nCols];
		//hack:  FIX!
		if (seismicVolume == null)
			seismicVolume = (StsSeismicVolume) currentModel.getCurrentObject(StsSeismicVolume.class);
		float absSeismicDataMax = Math.min(Math.abs(seismicVolume.dataMin), Math.abs(seismicVolume.dataMax));
		minDataAmplitude = minAmpFraction * absSeismicDataMax;
		initializeParameters();
		initializeToBoundingBox(croppedBoundingBox);
		initializeSliceInterpolation();

		int croppedRowMin = croppedBoundingBox.rowMin;
		int croppedRowMax = croppedBoundingBox.rowMax;
		int croppedColMin = croppedBoundingBox.colMin;
		int croppedColMax = croppedBoundingBox.colMax;
		int croppedSliceMin = croppedBoundingBox.sliceMin;
		int croppedSliceMax = croppedBoundingBox.sliceMax;
		int nCroppedSlices = croppedSliceMax - croppedSliceMin + 1;
		int nVolSlices = seismicVolume.nSlices;
		float[] traceValues = new float[nCroppedSlices];

		float[] nullTrace = new float[nCroppedSlices];
		Arrays.fill(nullTrace, seismicVolume.userNull);

		try
		{
			// row & col refer to the row and col in a croppedVolume over which picker is to run
			// volRow & volCol define the actual row and col in the volume (used only for reference)
			for (row = 0, volRow = croppedRowMin; volRow <= croppedRowMax; row++, volRow++)
			{
				//statusArea.setProgress(row*40.f/nRows);
				TracePoints[] prevRowTracesPoints = rowTracePoints;
				rowTracePoints = new TracePoints[nCols];
				TracePoints prevRowTracePoints = null;
				incrementLists();
				// get row plane from seismicVolume at this volRow
				FloatBuffer rowFloatBuffer = seismicVolume.getRowPlaneFloatBuffer(volRow, croppedColMin);
				if (rowFloatBuffer == null) return;
				// if(croppedColMin > 0) rowFloatBuffer.position(croppedColMin * nVolSlices);
				for (col = 0, volCol = croppedColMin; volCol <= croppedColMax; col++, volCol++)
				{
					// StsException.systemDebug(this, "constructPatchVolume", "col loop, col: " + col);
					rowFloatBuffer.position(volCol * nVolSlices + croppedSliceMin);
					rowFloatBuffer.get(traceValues);

					TracePoints tracePoints = null;
					if(!Arrays.equals(traceValues, nullTrace))
						tracePoints = new TracePoints(this, row, col, traceValues);
					rowTracePoints[col] = tracePoints;
					if(tracePoints == null) continue;
						// prevColTracePoints are tracePoints in prev row & same col
					TracePoints prevColTracePoints = null;
					if (prevRowTracesPoints != null)
						prevColTracePoints = prevRowTracesPoints[col];

					// here we add the connected patchPoints
					tracePoints.connectWindows(prevColTracePoints, prevRowTracePoints);

					prevRowTracePoints = tracePoints;
				}
				processPrevRowGrids(row);
				if (progressPanel == null) continue;
				if (progressPanel.isCanceled())
				{
					progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
					return;
				}
				progressPanel.setValue(row + 1);
			}
			addRemainingGrids();
			finish();
			getPatchVolumeClass().setDisplayCurvature(false);
			progressPanel.finished();
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "constructPatchVolume", e);
		}
	}

	public boolean setFilter()
	{
		return filter;
	}

	public void setFilter(boolean filter)
	{
		this.filter = filter;
	}

	public int getBoxFilterWidth()
	{
		return boxFilterWidth;
	}

	public void setBoxFilterWidth(int boxFilterWidth)
	{
		this.boxFilterWidth = boxFilterWidth;
	}

	protected StsPatchGrid mergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint)
	{
		StsPatchGrid mergedGrid, removedGrid;

		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

		if (otherPatchGrid.id < newPatchGrid.id)
		{
			mergedGrid = otherPatchGrid;
			removedGrid = newPatchGrid;
		}
		else
		{
			mergedGrid = newPatchGrid;
			removedGrid = otherPatchGrid;
		}
		// merge shouldn't fail, so a return of false indicates a problem: bail out
		if (!mergedGrid.mergePatchPoints(removedGrid))
		{
			StsException.systemError(this, "mergePatchGrids", "Failed to merge removedGrid " + removedGrid.toString() + " to mergedGrid " + mergedGrid.toString());
			return null;
		}

		//checkAddPatchGridToRowGrids(mergedGrid);
		//mergedGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
		removePatchGridFromLists(removedGrid);
		return mergedGrid;
	}
/*
	StsPatchGrid cantMergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
		// can't merge: find larger grid and add a clone of the overlapped window from the smaller grid to it

		StsPatchGrid changedGrid;
		PatchPoint clonedPoint;
		if (otherPatchGrid.nPatchPoints >= newPatchGrid.nPatchPoints)
		{
			changedGrid = otherPatchGrid;
			clonedPoint = newPatchPoint.clone();
			changedGrid.addPatchPoint(clonedPoint);
			changedGrid.addCorrelation(otherPatchPoint, clonedPoint, correl);
		}
		else
		{
			changedGrid = newPatchGrid;
			clonedPoint = otherPatchPoint.clone();
			changedGrid.addPatchPoint(clonedPoint);
			changedGrid.addCorrelation(clonedPoint, newPatchPoint, correl);
		}
		checkAddPatchGridToRowGrids(changedGrid);
		return changedGrid; // return null indicating this grid has been completely processed
	}
*/

	protected void checkAddPatchGridToRowGrids(StsPatchGrid patchGrid)
	{
		int patchID = patchGrid.id;
		if (patchGrid.rowGridAdded) return;
		StsPatchGrid value = rowGrids.put(patchID, patchGrid); // if return is null, no value exists at this key
		patchGrid.rowGridAdded = true;
		if (debug && StsPatchGrid.debugPatchID == patchGrid.id)
		{
			if (value == null)
				StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " added to rowGrids for row: " + row + " col: " + col);
			else
				StsException.systemDebug(this, "checkAddPatchGridToRowGrids", "patch " + patchID + " already exists for row: " + row + " col: " + col);
		}
	}

	private void removePatchGridFromLists(StsPatchGrid patchGrid)
	{
		StsPatchGrid value;
		int patchID = patchGrid.id;
		boolean debug = StsPatchGrid.debugPatchGrid && patchID == StsPatchGrid.debugPatchID;
		value = prevRowGrids.remove(patchID);
		if (debug)
		{
			if (value != null)
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from prevRowGrids for row: " + row);
			else
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in prevRowGrids for row: " + row);
		}
		value = rowGrids.remove(patchID);
		if (debug)
		{
			if (value != null)
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " removed from rowGrids for row: " + row);
			else
				StsException.systemDebug(this, "removePatchGridInGridList", "patch " + patchID + " doesn't exist in rowGrids for row: " + row);
		}
	}

	public boolean initialize(StsModel model)
	{
		initializeColorscale();
		initializeSliceInterpolation();
		initializePatchPointTotal();
		return true;
	}

	private void initializeSliceInterpolation()
	{
		nInterpolationIntervals = StsTraceUtilities.computeInterpolationInterval(zInc, 5);
		interpolatedZInc = zInc / nInterpolationIntervals;
		interpolatedSliceMin = 0;
		interpolatedSliceMax = (nSlices - 1) * nInterpolationIntervals;
		nInterpolatedSlices = interpolatedSliceMax + 1;
	}

	public void initialize()
	{
		clearSelectedPatches();
	}

	public void initializeColorscale()
	{
		try
		{
			if (colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				colorscale = new StsColorscale("Curvature", spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW), dataMin, dataMax);
				colorscale.setEditRange(dataMin, dataMax);
				colorscale.addActionListener(this);
			}
			colorscale.setRange(dataMin, dataMax);
		}
		catch (Exception e)
		{
			StsException.outputException("StsPatchcVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

	public void actionPerformed(ActionEvent e)
	{
		if (e.getSource() instanceof StsColorscale)
		{
			int modifiers = e.getModifiers();
			currentModel.viewObjectChangedAndRepaint(this, this);
		}
		else
		{
			//fireActionPerformed(e);
			currentModel.viewObjectChangedAndRepaint(this, this);
		}
		return;
	}

	private void initializeLists()
	{
		gridList = new ArrayList<>(100);
		rowGrids = new HashMap<>();
		rowGridsIterator = rowGrids.values().iterator();
		prevRowGrids = new HashMap<>(); // not used on first row
	}

	private void incrementLists()
	{
		prevRowGrids = rowGrids;
		prevRowGridsIterator = rowGridsIterator;
		clearPrevRowGridsAddedFlags();
		rowGrids = new HashMap<>();
		rowGridsIterator = rowGrids.values().iterator();
	}

	private void initializeParameters()
	{
		//maxStretchLimit = 1 + maxStretch;
		//minStretchLimit = 1/(1 + maxStretch);

		if (windowSize == 0) // windowSize is 0.5 wavelengths, half-window is 0.25 wavelengths
		{
			nHalfSamples = 1;
			windowEndIsZeroCrossing = true;
			halfWindowSize = 0.25f;
			halfWindowPickDifFactor = pickDifWavelengths / 0.25f;
		}
		else
		{
			boolean isHalfWave = !StsMath.isEven(windowSize);
			if (isHalfWave)
			{
				// window size is odd, so window ends with zero-crossing; half-window size is windowSize/2.
				// we need to find (windowSize +1)/2 zero-crossings above and below window center (which is a max or min).
				nHalfSamples = (windowSize + 1) / 2;
				windowEndIsZeroCrossing = true;
				halfWindowPickDifFactor = pickDifWavelengths * 2 / windowSize;
			}
			else
			{
				// window size is even, so window ends with same window type as center; half-window size is windowSize/2.
				// we need to find windowSize/2 points above and below with same window type as window center (which is a max or min).
				nHalfSamples = windowSize / 2;
				halfWindowPickDifFactor = pickDifWavelengths / nHalfSamples;
				windowEndIsZeroCrossing = false;
			}
		}
	}

	public StsPatchVolumeClass getPatchVolumeClass()
	{
		if (patchVolumeClass != null) return patchVolumeClass;
		patchVolumeClass = (StsPatchVolumeClass) getCreateStsClass();
		return patchVolumeClass;
	}

	private void finish()
	{
		// build patchGrid arrays. Overlapped grids will be split out into new grids
		// So construct an array of existing patchGrids and construct the various 2D arrays for each.
		// New grids generated will be added to the gridList
		int nGrids = gridList.size();
		rowSortedPatchGrids = new StsPatchGrid[nGrids];
		gridList.toArray(rowSortedPatchGrids);
		// debugCheckEmptyFraction(rowSortedPatchGrids);
		StsException.systemDebug(this, "finish", " Number of parent grids: " + nParentGrids + "Number of child grids: " + nGrids + " too small: " + nSmallGridsRemoved);
		// StsException.systemDebug(this, "finish", "max grid dimension: " + StsPatchGrid.maxGridSize);

		StsPatchGrid.sortRowFirst = true;
		Arrays.sort(rowSortedPatchGrids);
		// compute the total number of points
		initializePatchPointTotal();
		StsException.systemDebug(this, "finish", "Total number of points: " + nPointsTotal);
		// reset the index of each patch to its sequence in row-ordered array
		nGrids = 0;
		for (StsPatchGrid grid : rowSortedPatchGrids)
			grid.resetIndex(nGrids++);

		int num = getPatchVolumeClass().getSize();
		String newname = seismicName + ".patchVolume" + (num + 1);
		setName(newname);
		clearConstructionArrays();
		if (!isPersistent())
		{
			currentModel.getDatabase().blockTransactions();
			addToModel();
			currentModel.getDatabase().saveTransactions();
		}
		getPatchVolumeClass().setIsVisible(true);
		setIsVisible(true);
	}

	/*
		private void debugCheckEmptyFraction(StsPatchGrid[] grids)
		{
			float nTotal = 0;
			int[] used = null;
			float nUsed = 0;
			float nActualUsed = 0;

			for(StsPatchGrid grid : grids)
			{
				nTotal += grid.getGridSize();
				used  = grid.getGridPointsUsed();
				nUsed += used[0];
				nActualUsed += used[1];
			}
			float fractionUsed = nUsed/nTotal;
			float fractionActualUsed = nActualUsed/nTotal;
			StsException.systemDebug(this, "debugCheckEmptyFraction", "Fraction used: " + fractionUsed + "Fraction actualUsed: " + fractionActualUsed);
		}
	*/
	private void clearConstructionArrays()
	{
		gridList = null;
		rowGrids = null;
		prevRowGrids = null;
	}

	private void initializePatchPointTotal()
	{
		nPointsTotal = 0;
		if (rowSortedPatchGrids == null) return;
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			nPointsTotal += patchGrid.nPatchPoints;
	}

	private boolean checkColSortedPatchGrids()
	{
		if (colSortedPatchGrids != null) return true;
		if (rowSortedPatchGrids == null) return false;
		int nGrids = rowSortedPatchGrids.length;
		if (nGrids == 0) return false;
		colSortedPatchGrids = new StsPatchGrid[nGrids];
		System.arraycopy(rowSortedPatchGrids, 0, colSortedPatchGrids, 0, nGrids);
		StsPatchGrid.sortRowFirst = false;
		Arrays.sort(colSortedPatchGrids);
		return true;
	}

	private boolean checkRowSortedPatchGrids()
	{
		return rowSortedPatchGrids != null;
	}

	StsPatchGrid getGrid(int target)
	{
		int number = gridList.size();
		int high = number, low = -1, probe;
		while (high - low > 1)
		{
			probe = (high + low) / 2;
			int id = gridList.get(probe).idFinal;
			if (id > target)
				high = probe;
			else
				low = probe;
		}
		if (low == -1 || gridList.get(low).idFinal != target)
			return null;
		else
			return gridList.get(low);
	}

	private void clearPrevRowGridsAddedFlags()
	{
		Iterator<StsPatchGrid> prevRowGridIterator = prevRowGrids.values().iterator();
		while (prevRowGridIterator.hasNext())
		{
			StsPatchGrid patchGrid = prevRowGridIterator.next();
			patchGrid.rowGridAdded = false;
		}
	}

	/** Called for the last row only; unconditionally add all rowGrids unless they are too small. */
	void addRemainingGrids()
	{
		StsPatchGrid[] patchGrids = rowGrids.values().toArray(new StsPatchGrid[0]);
		for (int n = 0; n < patchGrids.length; n++)
		{
			StsPatchGrid patchGrid = patchGrids[n];
			if (!patchGrid.isTooSmall(nPatchPointsMin))
			{
				patchGrid.finish();
				gridList.add(patchGrid);
			}
		}
		StsException.systemDebug(this, "addRemainingGrids", "added " + patchGrids.length + " grids remaining.");
	}

	public void setCroppedBoundingBox(StsCroppedBoundingBox croppedBoundingBox)
	{
		this.croppedBoundingBox = croppedBoundingBox;
		croppedBoundingBox.setCroppedBoxRange();
	}

	/*
	* run the curvature calculation on the patches
	*/
	public void runCurvature(StsProgressPanel progressPanel, int filterSize, byte curveType, boolean runAllPatches)
	{
		this.filterSize = filterSize;

		if (runTimer)
		{
			timer = new StsTimer();
			timer.start();
		}
		StsPatchGrid[] runPatches = getRunCurvaturePatches(runAllPatches);
		int numPatches = runPatches.length;
		if (progressPanel != null)
			progressPanel.initialize(numPatches);
		int nValuePoints = 0;
		int nPoints = 0;
		float[] values = new float[nPointsTotal];
		double sum = 0;
		int progressUpdateInterval = Math.max(numPatches / 200, 1);
		int minNPoints = Math.min(filterSize * filterSize, StsQuadraticCurvature.minNPoints);
		for (int i = 0; i < numPatches; i++)
		{
			StsPatchGrid patch = runPatches[i];
			if (patch == null) continue;
			nPoints += patch.nPatchPoints;
			if (patch.computeCurvature(xInc, yInc, curveType, filterSize, minNPoints))
			{
				for (int row = patch.rowMin; row <= patch.rowMax; row++)
				{
					for (int col = patch.colMin; col <= patch.colMax; col++)
					{
						int ptCol = col - patch.colMin;
						int ptRow = row - patch.rowMin;
						float value = patch.curvature[ptRow][ptCol];
						if (value == StsPatchVolume.nullValue) continue;
						if (value == badCurvature || value == -badCurvature) continue;
						values[nValuePoints++] = value;
						sum += value;
					}
				}
			}

			if (progressPanel == null) continue;
			if (progressPanel.isCanceled())
			{
				progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
				clearPatches();
				return;
			}
			if (i % progressUpdateInterval == 0) progressPanel.setValue(i);
		}

		values = (float[]) StsMath.trimArray(values, nValuePoints);

		StsMessageFiles.infoMessage("Number curvature points: " + nValuePoints + " number of patch points " + nPointsTotal);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "Number curvature points: " + nValuePoints + " number of patch points " + nPointsTotal);


		double mean = sum / nValuePoints;
		int histogramDataInc = Math.max(1, nValuePoints / nHistogramValues);
		nHistogramValues = nValuePoints / histogramDataInc;
		histogramValues = new float[nHistogramValues + 1];
		double avgDev = 0;
		double sumSqr = 0;
		nHistogramValues = 0;
		for (int n = 0; n < nValuePoints; n++)
		{
			float value = values[n];
			double dif = value - mean;
			avgDev += Math.abs(dif);
			sumSqr += dif * dif;
			if (n % histogramDataInc == 0)
				histogramValues[nHistogramValues++] = value;
		}
		avgDev = avgDev / nValuePoints;
		double variance = (sumSqr) / nValuePoints;
		double stdDev = Math.sqrt(variance);
		StsMessageFiles.infoMessage("mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "mean " + mean + " avg dev: " + avgDev + " std dev: " + stdDev);
		dataMin = (float) (mean - 2.0 * avgDev);
		dataMax = (float) (mean + 2.0 * avgDev);
		colorscale.setRange(dataMin, dataMax);
		StsMessageFiles.infoMessage("range set to +- 2.0*avg dev: " + dataMin + " to " + dataMax);
		if (debug)
			StsException.systemDebug(this, "runCurvature", "range set to +- 2.0*std dev: " + dataMin + " to " + dataMax);

//        // reset outliers to dataMin/dataMax
//        for (StsPatchGrid patch : rowSortedPatchGrids)
//        {
//            if (patch == null) continue;
//            for (int row = patch.rowMin; row <= patch.rowMax; row++)
//            {
//                for (int col = patch.colMin; col <= patch.colMax; col++)
//                {
//                    float value = patch.getCurvature(row, col, dataMin, dataMax);
//                    if (value == StsPatchVolume.nullValue) continue;
//                }
//            }
//        }
		calculateHistogram(histogramValues, nHistogramValues);
		progressPanel.finished();

		if (runTimer) timer.stopPrint("Time to compute curvature for " + numPatches + " patches.");
	}

	private void clearPatches()
	{
		for (StsPatchGrid patch : rowSortedPatchGrids)
		{
			if (patch != null) patch.clear();
		}
	}

	private StsPatchGrid[] getRunCurvaturePatches(boolean runAllPatches)
	{
		if (runAllPatches || this.selectedPatchGrids == null) return rowSortedPatchGrids;
		else return selectedPatchGrids;
	}

	private ArrayList<TracePoints> getOtherTraces(TracePoints newTrace, TracePoints prevColTrace, TracePoints[] prevRowTraces)
	{
		ArrayList<TracePoints> otherTraces = new ArrayList<>();
		if (prevColTrace != null)
		{
			addOtherTrace(otherTraces, prevColTrace);
		}
		int col = newTrace.col;
		if (prevRowTraces != null)
		{
			// if (col > colMin)
			//    addOtherTrace(otherTraces, prevRowTraces[col - 1]);
			addOtherTrace(otherTraces, prevRowTraces[col]);
			// if (col < colMax)
			//    addOtherTrace(otherTraces, prevRowTraces[col + 1]);
		}
		return otherTraces;
	}

	private void addOtherTrace(ArrayList<TracePoints> otherTraces, TracePoints otherTrace)
	{
		if (otherTrace.nTracePatchPoints == 0) return;
		otherTraces.add(otherTrace);
	}

	private StsPatchGrid getPatchGrid(int id)
	{
		StsPatchGrid patchGrid = rowGrids.get(id);
		// if patchGrid exists in rowGrids, then it has already been added there and deleted from prevRowGrids
		if (patchGrid != null)
		{
			if (StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
				StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
						" gotten from rowGrids at row: " + row + " col: " + col);
			return patchGrid;
		}
		// if patchGrid is not in rowGrids, then add it there and delete it from prevRowGrids
		else if (prevRowGrids != null)
		{
			patchGrid = prevRowGrids.get(id);
			if (patchGrid != null)
			{
				rowGrids.put(id, patchGrid);
				prevRowGrids.remove(id);
				if (StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID))
					StsException.systemDebug(this, "getPatchGrid", "patch grid " + id +
							" gotten and deleted from prevRowsGrids and added to rowGrids at row: " + row + " col: " + col);
				return patchGrid;
			}
		}
		StsException.systemError(this, "getPatchGrid", "Couldn't get patchGrid for id " + id + " at row: " + " col: " + col);
		return null;
	}

	/**
	 * This PatchGridSet is for the row before row just finished.
	 * If a grid in this prev row is disconnected (doesn't have same patch in row just finished),
	 * then either delete it if it is a small window, or add it to volume set.
	 */
	void processPrevRowGrids(int row)
	{
		if (row == 0) return;
		StsPatchGrid[] prevRowPatchGrids = prevRowGrids.values().toArray(new StsPatchGrid[0]);
		int nDisconnectedGrids = 0;
		for (StsPatchGrid patchGrid : prevRowPatchGrids)
		{
			boolean disconnected = patchGrid.isDisconnected(row);
			if (!disconnected) continue;
			if (patchGrid.isTooSmall(nPatchPointsMin))
				nSmallGridsRemoved++;
			else
			{
				patchGrid.finish();
				gridList.add(patchGrid);
				nDisconnectedGrids++;
			}
			prevRowGrids.remove(patchGrid.id);
		}
		StsException.systemDebug(this, "processPrevRowGrids", "prev row: " + (row - 1) + " added " + nDisconnectedGrids + " disconnected grids");
	}

	public StsColorscale getCurvatureColorscale()
	{
		return colorscale;
	}

	/* Draw any map edges on all 2d sections */
	public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped,
							   boolean xAxisReversed, boolean yAxisReversed)
	{
		if (!getIsVisible()) return;

		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		StsColor drawColor = StsColor.BLACK; //getStsColor();
		gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
		drawColor.setGLColor(gl);

		if (dirNo == StsCursor3d.XDIR) /* map edge is along a col	*/
		{
			if (!checkColSortedPatchGrids()) return;
			int col = getNearestColCoor(dirCoordinate);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glShadeModel(GL.GL_SMOOTH);

			float x = dirCoordinate;
			// int nFirst = -1;
			// int n = -1;

			for (StsPatchGrid patchGrid : colSortedPatchGrids)
			{
				if (patchGrid.colMin > col) break;
				// n++;
				if (patchGrid.colMax < col) continue;

				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawCol(gl, col, x, yMin, yInc, colorscale, false, displayCurvature);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				// if (nFirst == -1) nFirst = n;
			}
			gl.glEnable(GL.GL_LIGHTING);
		}
		else if (dirNo == StsCursor3d.YDIR)
		{
			if (!checkRowSortedPatchGrids()) return;
			int row = getNearestRowCoor(dirCoordinate);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			gl.glShadeModel(GL.GL_SMOOTH);

			float y = dirCoordinate;
			// int nFirst = -1;
			// int n = -1;
			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			{
				if (patchGrid.rowMin > row) break;
				// n++;
				if (patchGrid.rowMax < row) continue;
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawRow(gl, row, y, xMin, xInc, colorscale, false, displayCurvature, filter, boxFilterWidth);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				// if (nFirst == -1) nFirst = n;
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) break;
			}
			gl.glEnable(GL.GL_LIGHTING);
		}
	}

	/** Draw any map edges on section */
	public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate)
	{
		if (!getIsVisible()) return;
		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		try
		{
			if (dirNo == StsCursor3d.ZDIR)
			{
				if (getDisplaySurfs())
					// displayPatchesNearXYZCursors(glPanel3d);
					displayPatchesNearZCursor(glPanel3d, dirCoordinate);
				return;
			}
			if (dirNo == StsCursor3d.YDIR)
			{
				if (!checkRowSortedPatchGrids()) return;
				gl.glDisable(GL.GL_LIGHTING);
				gl.glShadeModel(GL.GL_SMOOTH);
				StsColor drawColor = StsColor.BLACK; //getStsColor();
				drawColor.setGLColor(gl);
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
				int row = getNearestRowCoor(dirCoordinate);
				if (row == -1) return;
				float xMin = getXMin();
				float xInc = getXInc();
				for (StsPatchGrid patchGrid : rowSortedPatchGrids)
				{
					if (patchGrid.rowMin > row) break;
					if (patchGrid.rowMax < row) continue;
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawRow(gl, row, dirCoordinate, xMin, xInc, colorscale, true, displayCurvature, setFilter(), getBoxFilterWidth());
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
					gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				}
			}
			else if (dirNo == StsCursor3d.XDIR)
			{
				if (!checkColSortedPatchGrids()) return;
				gl.glDisable(GL.GL_LIGHTING);
				gl.glShadeModel(GL.GL_SMOOTH);
				StsColor drawColor = StsColor.BLACK; //getStsColor();
				drawColor.setGLColor(gl);
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
				int col = getNearestColCoor(dirCoordinate);
				if (col == -1) return;
				float yMin = getYMin();
				float yInc = getYInc();
				for (StsPatchGrid patchGrid : colSortedPatchGrids)
				{
					if (patchGrid.colMin > col) break;
					if (patchGrid.colMax < col) continue;
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawCol(gl, col, dirCoordinate, yMin, yInc, colorscale, true, displayCurvature);
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchID) ;
					gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				}
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "drawOnCursor3d", e);
		}
		finally
		{
			glPanel3d.resetViewShift(gl);
			gl.glEnable(GL.GL_LIGHTING);
		}

	}

	public void pickOnCursor3d(StsGLPanel3d glPanel3d)
	{
		StsCursor3d cursor3d = glPanel3d.window.getCursor3d();
		if (cursor3d == null) return;
		GL gl = glPanel3d.getGL();
		for (int dir = 0; dir < 2; dir++)
		{
			float dirCoordinate = cursor3d.getCurrentDirCoordinate(dir);
			drawOnCursor3d(glPanel3d, dir, dirCoordinate);
		}
	}

	public void display(StsGLPanel glPanel)
	{
		if (!getDisplaySurfs() || selectedPatchGrids == null) return;
		GL gl = glPanel.getGL();

		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();

		initializePatchDraw(gl);
		for (StsPatchGrid patchGrid : selectedPatchGrids)
			patchGrid.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);

		if (getDisplayVoxels())
		{
			displayVoxels(glPanel);
		}
	}

	public void displayVoxelsCursor(StsGLPanel glPanel3d, StsPoint[] points, boolean is3d)
	{
		//System.out.println("Display Voxels");
		GL gl = glPanel3d.getGL();
		if (gl == null) return;
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		StsGridPoint point1 = new StsGridPoint(points[0], this);
		StsGridPoint point2 = new StsGridPoint(points[3], this);
		int sameRow = StsGridPoint.getSameRow(point1, point2); // if not -1, this is the row these two points are on
		if (sameRow != -1)
		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
			gl.glColor4f(1.f, 1.f, 1.f, 1.f);

			float y;

			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			{
				int row1 = sameRow - 4;
				int row2 = sameRow + 4;
				y = yMin + (yInc * row1);
				for (int n = row1; n <= row2; n++)
				{
					//if (n == 146)
					//System.out.println("draw vox row "+n+" "+patchGrid.colMin+" "+patchGrid.colMax+" "+y);
					patchGrid.drawRow(gl, n, y, xMin, xInc, colorscale, is3d, displayCurvature, setFilter(), getBoxFilterWidth());
					y += yInc;
				}
			}

			glPanel3d.resetViewShift(gl);
			gl.glEnable(GL.GL_LIGHTING);
			return;
		}

	}

	public void displayVoxels(StsGLPanel glPanel3d)
	{
		//System.out.println("Display Voxels");

		GL gl = glPanel3d.getGL();
		if (gl == null) return;

		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
			gl.glColor4f(1.f, 1.f, 1.f, 1.f);
			float xMin = getXMin();
			float xInc = getXInc();
			float yMin = getYMin();
			float yInc = getYInc();
			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
				patchGrid.drawRowVox(gl, yMin, yInc, xMin, xInc, colorscale);
			gl.glEnable(GL.GL_LIGHTING);
			return;
		}

	}

	// surfaces near cursor only to keep clutter down
	public void displayPatchesNearZCursor(StsGLPanel glPanel3d, float z)
	{
		//System.out.println("Display Surfaces");
		GL gl = glPanel3d.getGL();
		if (gl == null) return;

		initializePatchDraw(gl);
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
		{
			if (patchGrid.isPatchGridNearZCursor(z))
				patchGrid.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);
		}
		return;
	}

	private void initializePatchDraw(GL gl)
	{
		gl.glEnable(GL.GL_LIGHTING);
		gl.glShadeModel(GL.GL_SMOOTH);
		gl.glEnable(GL.GL_NORMALIZE);
		gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE);
		// gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL, GL.GL_SEPARATE_SPECULAR_COLOR);
		//gl.glDisable(GL.GL_CULL_FACE);
		//gl.glDisable(GL.GL_POLYGON_STIPPLE);
		gl.glColor4f(1.f, 1.f, 1.f, 1.f);
	}

	public void displayPatchesNearXYZCursors(StsGLPanel glPanel3d)
	{
		StsCursor3d cursor3d = currentModel.win3d.getCursor3d();
		GL gl = glPanel3d.getGL();
		if (gl == null) return;

		float x = cursor3d.getCurrentDirCoordinate(XDIR);
		float y = cursor3d.getCurrentDirCoordinate(YDIR);
		float z = cursor3d.getCurrentDirCoordinate(ZDIR);
		StsPoint cursorPoint = new StsPoint(x, y, z);
		boolean displayCurvature = getPatchVolumeClass().getDisplayCurvature();
		if (currentCursorPoint != null && currentCursorPoint.equals(cursorPoint) && cursorPointPatch != null)
		{
			drawPatch(cursorPointPatch, displayCurvature, gl);
			return;
		}

		cursorPoint = null;
		currentCursorPoint = null;

		int volumeRow = getNearestRowCoor(y);
		if (volumeRow == -1) return;
		int volumeCol = getNearestColCoor(x);
		if (volumeCol == -1) return;
		int slice = getNearestSliceCoor(z);
		if (slice == -1) return;
		float dzPatch = largeFloat;
		cursorPointPatch = null;
		for (StsPatchGrid patchGrid : rowSortedPatchGrids)
		{
			if (patchGrid == null) continue;
			// if(patchGrid.values == null) continue;
			if (patchGrid.rowMin > volumeRow) break;
			if (patchGrid.rowMax >= volumeRow)
			{
				if (patchGrid.colMin <= volumeCol && patchGrid.colMax >= volumeCol)
				{
					float dz = patchGrid.getZDistance(volumeRow, volumeCol, z);
					if (dz < dzPatch)
					{
						cursorPointPatch = patchGrid;
						dzPatch = dz;
						currentCursorPoint = new StsPoint(x, y, z);
					}
				}
			}
		}
		if (cursorPointPatch == null) return;

		drawPatch(cursorPointPatch, displayCurvature, gl);

		return;
	}

	private void drawPatch(StsPatchGrid patch, boolean displayCurvature, GL gl)
	{
		initializePatchDraw(gl);
		patch.draw(gl, xMin, xInc, yMin, yInc, displayCurvature, colorscale);
	}

//	public void setWizard(StsVolumeCurvatureWizard wizard) {
//		this.wizard = wizard;
//	}
//
//	public StsVolumeCurvatureWizard getWizard() {
//		return wizard;
//	}

	public int[] getPatchRangeForRow(int volumeRow)
	{
		int rowMin = -1;
		int rowMax = -1;
		int nPatchGrids = rowSortedPatchGrids.length;
		int n = 0;
		for (n = 0; n < nPatchGrids; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			if (patchGrid == null) continue;
			rowMax = n - 1;
			if (rowMin == -1)
				if (patchGrid.rowMin <= volumeRow && patchGrid.rowMax >= volumeRow) rowMin = n;
				else if (patchGrid.rowMin > volumeRow)
					break;
		}
		if (rowMin == -1)
			return new int[]{0, 0};
		else
			return new int[]{rowMin, rowMax};
	}

	public int[] getPatchRangeForCol(int col)
	{
		int colMin = -1;
		int colMax = -1;
		int nPatchGrids = colSortedPatchGrids.length;
		for (int n = 0; n < nPatchGrids; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			if (patchGrid == null) continue;
			if (colMin == -1)
			{
				if (patchGrid.colMin <= col && patchGrid.colMax >= col) colMin = n;
			}
			else if (patchGrid.colMin > col)
			{
				colMax = n - 1;
				break;
			}
		}
		if (colMin == -1 || colMax == -1)
			return new int[]{0, 0};
		else
			return new int[]{colMin, colMax};
	}

	public boolean getTraceCurvature(int volRow, int volCol, float[] buffer, int[] patchRange)
	{
		try
		{
			Arrays.fill(buffer, nullValue);
			int nPatchMin = patchRange[0];
			int nPatchMax = patchRange[1];
			boolean traceLoaded = false;
			for (int nPatch = nPatchMin; nPatch <= nPatchMax; nPatch++)
			{
				StsPatchGrid patchGrid = rowSortedPatchGrids[nPatch];
				if (patchGrid == null) continue;
				if (patchGrid.curvature == null) continue;
				int patchRow = volRow - patchGrid.rowMin;
				int patchCol = volCol - patchGrid.colMin;
				if (patchRow < 0 || patchCol < 0 || patchRow >= patchGrid.nRows || patchCol >= patchGrid.nCols)
					continue;
				float[][] pointsZ = patchGrid.getPointsZ();
				if (pointsZ == null) continue;
				float z = pointsZ[patchRow][patchCol];
				if (z == StsParameters.nullValue) continue;
				float val = patchGrid.curvature[patchRow][patchCol];
				if (val == nullValue) continue;
				int slice = getNearestSliceCoor(z);
				buffer[slice] = val;
				traceLoaded = true;
			}
			return traceLoaded;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "getTraceCurvature", e);
			return false;
		}
	}

	public int getNearestSliceCoor(float z)
	{
		int slice = Math.round((z - zMin) / interpolatedZInc);
		if (slice < 0 || slice >= nInterpolatedSlices) return -1;
		return slice;
	}

	public int getPatchPointIndex(PatchPoint patchPoint)
	{
		return patchPoint.getRow() * nCols + patchPoint.getCol();
	}

	public String toString()
	{
		return name;
	}

	public void treeObjectSelected()
	{
		getPatchVolumeClass().selected(this);
		currentModel.getGlPanel3d().checkAddView(StsView3d.class);
		currentModel.win3dDisplayAll();
	}


	public Object[] getChildren()
	{
		return new Object[0];
	}

	public StsFieldBean[] getDisplayFields()
	{
		//displayAttributeBean.setListItems(displayAttributes);
		return displayFields;
	}

	public StsFieldBean[] getPropertyFields()
	{
		return propertyFields;
	}


	public StsObjectPanel getObjectPanel()
	{
		if (objectPanel == null)
		{
			objectPanel = StsObjectPanel.constructor(this, true);
		}
		return objectPanel;
	}

	public boolean anyDependencies()
	{
		return false;
	}

	public StsColorscale getColorscale()
	{
		//setDataHistogram();
		return colorscale;
	}

	public void setColorscale(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
		currentModel.win3dDisplayAll();
	}

	public void setDisplaySurfs(boolean displaySurfs)
	{
		if (this.displaySurfs == displaySurfs)
			return;
		this.displaySurfs = displaySurfs;
		currentModel.win3dDisplayAll();
	}

	public boolean getDisplaySurfs()
	{
		return displaySurfs;
	}

	public void setDisplayVoxels(boolean displayVoxels)
	{
		if (this.displayVoxels == displayVoxels)
			return;
		this.displayVoxels = displayVoxels;
		currentModel.win3dDisplayAll();
	}

	public boolean getDisplayVoxels()
	{
		return displayVoxels;
	}

	public void setIsVisible(boolean vis)
	{
		super.setIsVisible(vis);
		currentModel.win3dDisplayAll();
	}

	public void setDataMin(float min)
	{
		dataMin = min;
		if (colorscale == null) return;
		colorscale.setRange(dataMin, dataMax);
	}

	public void setDataMax(float max)
	{
		dataMax = max;
		if (colorscale == null) return;
		colorscale.setRange(dataMin, dataMax);
	}

	public void addRemoveSelectedPatch(StsCursorPoint cursorPoint)
	{
		boolean displayChildPatches = getPatchVolumeClass().getDisplayChildPatches();
		float[] xyz = cursorPoint.point.v;
		int volumeRow = getNearestRowCoor(xyz[1]);
		int volumeCol = getNearestColCoor(xyz[0]);
		float z = xyz[2];
		StsPatchGrid selectedPatch = getNearestPatch(volumeRow, volumeCol, z);
		if (selectedPatch == null) return;

		float iline = getRowNumFromRow(volumeRow);
		float xline = getColNumFromCol(volumeCol);
		StsMessageFiles.logMessage("Picked patch parent id: " + selectedPatch.id + " child id: " + selectedPatch.idFinal + " at iline: " + iline + " xline: " + xline);
	/*
		if(cursorPoint.dirNo == StsCursor3d.YDIR)
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeRowCorrel(volumeRow, volumeCol));
        else //dirNo == XDIR
            StsMessageFiles.logMessage("     volumeRow correl: " + selectedPatch.getVolumeColCorrel(volumeRow, volumeCol));
    */
		int parentID = selectedPatch.id;
		// int nSheet = selectedPatch.nSheet;
		String addRemove;
		boolean removePatch = StsMath.arrayContains(selectedPatchGrids, selectedPatch);
		if (removePatch)
			addRemove = " removed ";
		else
			addRemove = " added ";
		int nGridsAdded = 0;
		if (!displayChildPatches)
		{
			if (removePatch)
				selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayDeleteElement(selectedPatchGrids, selectedPatch);
			else
				selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayAddElement(selectedPatchGrids, selectedPatch);
		}
		else
		{
			for (StsPatchGrid patchGrid : rowSortedPatchGrids)
			{
				if (patchGrid.id == parentID) // && patchGrid.nSheet == nSheet)
				{
					nGridsAdded++;
					if (removePatch)
						selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayDeleteElement(selectedPatchGrids, patchGrid);
					else
						selectedPatchGrids = (StsPatchGrid[]) StsMath.arrayAddElement(selectedPatchGrids, patchGrid);
				}
			}
		}
		StsMessageFiles.logMessage("Picked patch parent id: " + selectedPatch.id + addRemove + nGridsAdded + " child grids.");
	}

	private void clearSelectedPatches()
	{
		selectedPatchGrids = null;
	}

	public StsPatchGrid getNearestPatch(int volumeRow, int volumeCol, float z)
	{
		StsPatchGrid nearestPatch = null;
		float nearestPatchZ = largeFloat;
		//TODO make an iterator which returns patches which cross this volumeRow
		int[] patchRange = this.getPatchRangeForRow(volumeRow);
		int patchMin = patchRange[0];
		int patchMax = patchRange[1];
		for (int n = patchMin; n <= patchMax; n++)
		{
			StsPatchGrid patchGrid = rowSortedPatchGrids[n];
			float patchZ = patchGrid.getPointZ(volumeRow, volumeCol);
			if (patchZ == nullValue) continue;
			float dz = Math.abs(z - patchZ);
			if (dz < nearestPatchZ)
			{
				nearestPatch = patchGrid;
				nearestPatchZ = dz;
			}
		}
		return nearestPatch;
	}
/*
    public void addSelectedPatch(int patchID)
    {
        if(patchID < 0 || patchID >= rowSortedPatchGrids.length)
        {
            StsException.systemError(this, "addSelectedPatch", "patchID out of range: " + patchID);
            return;
        }
        StsPatchGrid selectedPatch = this.rowSortedPatchGrids[patchID];
        selectedPatchGrids = (StsPatchGrid[])StsMath.arrayAddElement(selectedPatchGrids, selectedPatch);
    }
*/
}

class TracePoints
{
	StsPatchVolume patchVolume;
	/** volume row for this trace */
	int row;
	/** volume col for this trace */
	int col;
	/** array of PatchPoints for this trace of various pointTypes (min, max, +zero-crossing, -zero-crossing; can be false or missing) */
	PatchPoint[] tracePatchPoints = new PatchPoint[0];
	/** length of tracePatchPoints array */
	int nTracePatchPoints;
	/** offset to first plus-zero-crossing in tracePatchPoints array */
	int zeroPlusOffset;
	/** a half-wave length window for each legitimate window (acceptable pointType) */
	CorrelationWindow[] windows;
	/** number of windows in this trace; one for which tracePatchPoint */
	int nWindows;
	/** double-linked list of connections between this trace and prevCol */
	ConnectionList colConnections;
	/** double-linked list of connections between this trace and prevCol */
	ConnectionList rowConnections;
	/** closest connections between this trace and the prevColTrace */
	CorrelationWindow[] colCloseConnectWindows;
	/** closest connections between this trace and the prevRowTrace */
	CorrelationWindow[] rowCloseConnectWindows;
	/**
	 * From the traceValues array, create an array of PatchPoint for each traceValue which qualifies as a pickType (all, min, max, +/- zero-crossing.
	 * This will be a sequential subset of the values array with nTracePatchPoints in this tracePatchPoints array
	 * @param patchVolume
	 * @param row volume row of this trace
	 * @param col volume col of this trace
	 * @param traceValues original seismic values for this trace
	 */
	TracePoints(StsPatchVolume patchVolume, int row, int col, float[] traceValues)
	{
		this.patchVolume = patchVolume;
		this.row = row;
		this.col = col;
		// tracePoints - uniform cubic interpolation of traceValues
		float[] tracePoints = StsTraceUtilities.computeCubicInterpolatedPoints(traceValues, patchVolume.nInterpolationIntervals);
		float z = patchVolume.croppedBoundingBox.zMin;
		if (tracePoints == null) return;
		int nTracePoints = tracePoints.length;
		tracePatchPoints = new PatchPoint[nTracePoints];
		int nTracePatchPoint = 0;
		byte[] tracePointTypes = StsTraceUtilities.getPointTypes(tracePoints);

		// create tracePoints from values and pointTypes for legitimate events (zero+, max, zero-, min)
		// count the number missing so we have the final array length with missing points added
		int nTotalMissing = 0;
		int nMissing;
		PatchPoint prevPoint, nextPoint = null;

		try
		{
			for (int n = 0; n < nTracePoints; n++, z += patchVolume.interpolatedZInc)
			{
				byte tracePointType = tracePointTypes[n];
				if (StsTraceUtilities.isMaxMinZeroOrFalseMaxMin(tracePointType))
				{

					prevPoint = nextPoint;
					nextPoint = new PatchPoint(this, row, col, n, z, tracePoints[n], tracePointType, nTracePatchPoint);
					tracePatchPoints[nTracePatchPoint] = nextPoint;
					nTracePatchPoint++;
					nMissing = getNMissingPoints(prevPoint, nextPoint);
					nTotalMissing += nMissing;
				}
			}
			// insert any missing points so that we have a series of zero+, max, zero-, min points with any missing filled as ZPM, MXM, MNM, ZMM

			nTracePatchPoints = nTracePatchPoint;
			if (nTotalMissing > 0)
			{
				int nTotalPoints = nTracePatchPoints + nTotalMissing;
				PatchPoint[] addedPoints = new PatchPoint[nTotalPoints];
				nextPoint = tracePatchPoints[0];
				nTracePatchPoint = 0;
				for (int n = 1; n < nTracePatchPoints; n++)
				{
					prevPoint = nextPoint;
					nextPoint = tracePatchPoints[n];

					addedPoints[nTracePatchPoint] = prevPoint.resetIndex(nTracePatchPoint);
					nTracePatchPoint++;
					nMissing = getNMissingPoints(prevPoint, nextPoint);
					if( nMissing > 0)
					{
						byte pointTypeStart = StsTraceUtilities.pointTypesAfter[prevPoint.pointType];
						byte pointTypeEnd = StsTraceUtilities.pointTypesBefore[nextPoint.pointType];
						byte missingType = pointTypeStart;
						for (int i = 0; i < nMissing; i++, missingType++)
						{
							if (missingType > 4) missingType -= 4;
							addedPoints[nTracePatchPoint] = new PatchPoint(prevPoint, missingType, nTracePatchPoint);
							nTracePatchPoint++;
						}
					}
				}
				addedPoints[nTracePatchPoint] = nextPoint.resetIndex(nTracePatchPoint);
				nTracePatchPoint++;
				//if (nTracePatchPoint != nTotalPoints)
				//	StsException.systemError(this, "new TracePoints", " nTotalPoints " + nTotalPoints + " not equal to nTracePatchPoints " + nTracePatchPoint);
				tracePatchPoints = addedPoints;
				nTracePatchPoints = nTracePatchPoint;
			}
			else
				tracePatchPoints = (PatchPoint[]) StsMath.trimArray(tracePatchPoints, nTracePatchPoints);

			zeroPlusOffset = StsTraceUtilities.zeroPlusOffset[tracePatchPoints[0].pointType];
			constructTraceWindows();
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "TracePoints(vol, row, col, values)", e);
		}
	}

	static public boolean pointsNotInSequence(PatchPoint prevPoint, PatchPoint nextPoint)
	{
		if(prevPoint == null || nextPoint == null) return false;
		byte prevType = StsTraceUtilities.coercedPointTypes[prevPoint.pointType];
		byte nextType = StsTraceUtilities.coercedPointTypes[nextPoint.pointType];
		return StsTraceUtilities.pointTypesAfter[prevType] != nextType;
	}

	static int getNMissingPoints(PatchPoint prevPoint, PatchPoint nextPoint)
	{
		if(prevPoint == null || nextPoint == null) return 0;
		byte prevType = StsTraceUtilities.coercedPointTypes[prevPoint.pointType];
		byte nextType = StsTraceUtilities.coercedPointTypes[nextPoint.pointType];
		if(StsTraceUtilities.pointTypesAfter[prevType] == nextType) return 0;
		byte pointTypeStart = StsTraceUtilities.pointTypesAfter[prevPoint.pointType];
		byte pointTypeEnd = StsTraceUtilities.pointTypesBefore[nextPoint.pointType];
		return StsTraceUtilities.getNumPointTypesBetweenInclusive(pointTypeStart, pointTypeEnd);
	}

	/**
	 * creates correlated connections between this trace and traces at prev row & same col and prev col & same row
	 * @param prevColTrace prev trace in same col, prev row
	 * @param prevRowTrace prev trace in same row, prev col
	 */
	void connectWindows(TracePoints prevColTrace, TracePoints prevRowTrace)
	{
		if (prevRowTrace == null && prevColTrace == null) return;

		int nTracePatchPoints = tracePatchPoints.length;
		if (nTracePatchPoints == 0) return;

		// The ConnectionLists for this trace is initialized with top and bot inactive connections to prev row and col traces.
		// These inactive connections are used to limit the search process.
		// a connection is from the first window in this patchPointsList to the first window in the corresponding row or col trace
		if(prevColTrace != null)
		{
			colConnections = new ConnectionList(this, prevColTrace);
			prevColTrace.colConnections = colConnections;
		}
		if(prevRowTrace != null)
		{
			rowConnections = new ConnectionList(this, prevRowTrace);
			prevRowTrace.rowConnections = rowConnections;
		}
		// create initial guesses of connections from the windows on this trace to windows on prevCol and prevRow traces
		// guesses are the closest windows vertically on the other traces to the window on this trace
		if(prevColTrace != null)
			createClosestConnectWindow(this, prevColTrace, false);
		if(prevRowTrace != null)
			createClosestConnectWindow(this, prevRowTrace, true);
		// Iterate from maxCorrelation down to minCorrelation in nIterations steps.  At each iteration,
		// make connections from prevCol and prevRow traces to this trace */
		for (int iter = 0; iter < patchVolume.nIterations; iter++)
			connectWindows(patchVolume, prevColTrace, prevRowTrace, iter);
	}

	/** This prevents crossing connections.
	 *  Create lists of closest points of same type from trace to otherTrace and back.
	 *  If connections are identical, then retain in the trace->otherTrace list; if not null out.
	 * @param trace connections will be from this trace to otherTrace
	 * @param otherTrace connections list back will be compared with trace list
	 * @param isRow indicates trace and otherTrace are on the same row (other trace is prevCol)
	 */
	static void createClosestConnectWindow(TracePoints trace, TracePoints otherTrace, boolean isRow)
	{
		CorrelationWindow window, otherWindow, backWindow;
		int n;
		CorrelationWindow[] windows = trace.windows;
		int nWindows = trace.nWindows;
		CorrelationWindow[] connectWindows = trace.createClosestConnectWindows(otherTrace);
		CorrelationWindow[] otherConnectWindows = otherTrace.createClosestConnectWindows(trace);
		try
		{
			CorrelationWindow prevWindow = null;
			for (n = 0; n < nWindows; n++)
			{
				window = windows[n];
				otherWindow = connectWindows[n];

				if(otherWindow == null)
					continue;
				int otherIndex = otherWindow.windowIndex;
				backWindow = otherConnectWindows[otherIndex];
				if (backWindow != window)
				{
					connectWindows[n] = null;
					otherConnectWindows[otherIndex] = null;
				}
				else if(prevWindow != null && otherWindow.centerPoint.slice <= prevWindow.centerPoint.slice)
					connectWindows[n] = null;
				else
					prevWindow = otherWindow;
			}
			// check for crossing connections; remove connection if it crosses
			if (isRow)
			{
				trace.rowCloseConnectWindows = connectWindows;
				otherTrace.rowCloseConnectWindows = otherConnectWindows;
			}
			else
			{
				trace.colCloseConnectWindows = connectWindows;
				otherTrace.colCloseConnectWindows = otherConnectWindows;
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(TracePoints.class, "createClosestConnectWindow", e);
		}
	}

	static void checkCrossings(CorrelationWindow[] connectWindows, int nextIndex)
	{
		if(nextIndex == 0) return;

		int prevSlice = -1, slice, nextSlice;
		if(nextIndex > 1)
			prevSlice = connectWindows[nextIndex-2].centerPoint.slice;
		slice = connectWindows[nextIndex-1].centerPoint.slice;
		nextSlice = connectWindows[nextIndex].centerPoint.slice;
		if(prevSlice >=slice || nextSlice <= slice)
			connectWindows[nextIndex-1] = null;
	}

	CorrelationWindow[] createClosestConnectWindows(TracePoints otherTrace)
	{
		CorrelationWindow[] closeConnectWindows = new CorrelationWindow[nWindows];

		// assign closest otherWindow to this window regardless of pointType
		CorrelationWindow[] otherWindows = otherTrace.windows;
		int nOtherWindows = otherTrace.nWindows;
		CorrelationWindow otherWindowAbove = otherWindows[0];
		CorrelationWindow otherWindowBelow = otherWindows[1];
		int otherNextIndex = 2;
		for (int i = 0; i < nWindows; i++)
		{
			CorrelationWindow window = windows[i];
			if (window.isBelowOrEqual(otherWindowAbove) && window.isAboveOrEqual(otherWindowBelow))
			{
				CorrelationWindow otherWindow = window.getClosestWindow(otherWindowAbove, otherWindowBelow);
				closeConnectWindows[i] = otherWindow;
			}
			else if (window.isAboveOrEqual(otherWindowAbove)) continue;
			else // window is below otherWindowBelow, so move otherWindows down
			{
				while (window.isBelowOrEqual(otherWindowBelow) && otherNextIndex < nOtherWindows)
				{
					otherWindowAbove = otherWindowBelow;
					otherWindowBelow = otherWindows[otherNextIndex++];
				}
				if (window.isBelowOrEqual(otherWindowAbove) && window.isAboveOrEqual(otherWindowBelow))
				{
					CorrelationWindow otherWindow = window.getClosestWindow(otherWindowAbove, otherWindowBelow);
					closeConnectWindows[i] = otherWindow;
				}
			}
		}
		// now adjust guesses to nearest otherWindow of this pointType
		for (int i = 0; i < nWindows; i++)
		{
			CorrelationWindow window = windows[i];
			CorrelationWindow otherWindow = closeConnectWindows[i];
			if(otherWindow == null) continue;
			int pointTypeOffset = getPointTypeDif(window, otherWindow);
			if(pointTypeOffset == 0) continue; // we have the correct type: so don't adjust
			if (pointTypeOffset >= 2) pointTypeOffset -= 4;
			int otherWindowIndex = otherWindow.windowIndex + pointTypeOffset;
			if (otherWindowIndex < 0)
				otherWindowIndex += 4;
			else if(otherWindowIndex >= nOtherWindows)
				otherWindowIndex -= 4;
			closeConnectWindows[i] = otherWindows[otherWindowIndex];
		}
		return closeConnectWindows;
	}

	private void connectWindows(StsPatchVolume patchVolume, TracePoints prevColTrace, TracePoints prevRowTrace, int iter)
	{
		CorrelationWindow matchingWindow, backMatchingWindow;
		Connection colConnection, rowConnection;
		CorrelationWindow otherClosestWindow, closestWindow;
		CorrelationWindow connectedWindowAbove, connectedWindowBelow;
		CorrelationWindow window;

		try
		{
			reinitializeTraceIndices(prevRowTrace, prevColTrace);
			for (int n = 0; n < nWindows; n++)
			{
				window = windows[n];
				colConnection = null;
				if (prevColTrace != null && !window.hasColConnection())
				{
					otherClosestWindow = colCloseConnectWindows[n];
					connectedWindowAbove = colConnections.connectionAbove.otherWindow;
					connectedWindowBelow = colConnections.connectionBelow.otherWindow;
					matchingWindow = TracePoints.connectWindows(patchVolume, window, otherClosestWindow, prevColTrace, connectedWindowAbove, connectedWindowBelow, iter);
					if (matchingWindow != null && patchVolume.backMatch)
					{
						closestWindow = prevColTrace.colCloseConnectWindows[matchingWindow.windowIndex];
						connectedWindowAbove = prevColTrace.colConnections.connectionAbove.window;
						connectedWindowBelow = prevColTrace.colConnections.connectionBelow.window;
						backMatchingWindow = TracePoints.connectWindows(patchVolume, matchingWindow, closestWindow, this, connectedWindowAbove, connectedWindowBelow, iter);
						if (backMatchingWindow != null && backMatchingWindow != window && backMatchingWindow.stretchCorrelation >= matchingWindow.stretchCorrelation)
							matchingWindow = null;
					}
					if (matchingWindow != null)
						colConnection = new Connection(matchingWindow, window);
				}
				rowConnection = null;
				if (prevRowTrace != null && !window.hasRowConnection())
				{
					otherClosestWindow = rowCloseConnectWindows[n];
					connectedWindowAbove = rowConnections.connectionAbove.otherWindow;
					connectedWindowBelow = rowConnections.connectionBelow.otherWindow;
					matchingWindow = TracePoints.connectWindows(patchVolume, window, otherClosestWindow, prevRowTrace, connectedWindowAbove, connectedWindowBelow, iter);
					if (matchingWindow != null && patchVolume.backMatch)
					{
						closestWindow = prevRowTrace.rowCloseConnectWindows[matchingWindow.windowIndex];
						connectedWindowAbove = prevRowTrace.rowConnections.connectionAbove.window;
						connectedWindowBelow = prevRowTrace.rowConnections.connectionBelow.window;
						backMatchingWindow = TracePoints.connectWindows(patchVolume, matchingWindow, closestWindow, this, connectedWindowAbove, connectedWindowBelow, iter);
						if (backMatchingWindow != null && backMatchingWindow != window && backMatchingWindow.stretchCorrelation >= matchingWindow.stretchCorrelation)
							matchingWindow = null;
					}
					if (matchingWindow != null)
						rowConnection = new Connection(matchingWindow, window);
				}
				processNewConnections(window, colConnection, rowConnection);
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "connectWindows(vol,trace,trace,iter)", e);
		}
	}

	/** from the current connectionAbove for this connectionList, find index offset to this window and apply it to the otherTrace window
	 *  to get the middle matchingWindow.  Check connections to otherTrace windows above and below.
	 *
	 * @param patchVolume
	 * @param newWindow newWindow we want to match on this otherTrace
	 * @param centerOtherWindow center candidate window for matching
	 * @param otherTrace other trace to which we want to make a connection from this trace
	 * @param otherConnectWindowAbove otherTrace windows on connections above
	 * @param otherConnectWindowBelow otherTrace windows on connections below
	 * @param iter iteration we are on
	 * @return best connection between newWindow and a matchingWindow on this trace; return null of non exists or don't qualify
	 */
	static private CorrelationWindow connectWindows(StsPatchVolume patchVolume, CorrelationWindow newWindow,
													CorrelationWindow centerOtherWindow, TracePoints otherTrace,
													CorrelationWindow otherConnectWindowAbove, CorrelationWindow otherConnectWindowBelow,
													int iter)
	{
		// matchingWindow we wish to find
		CorrelationWindow matchingWindow = null;
		// index of candidate centerOtherWindow
		int centerOtherWindowIndex;
		// indexes for the two other matchingWindow candidates above and below the center (offsets of -4 and +4)
		int aboveWindowIndex, belowWindowIndex;
		// candidate matching windows above and below
		CorrelationWindow aboveOtherWindow, belowOtherWindow;
		// having found a matchingWindow, try matching it back to the newTrace to see if we find a different connection with a better correlation
		// if we do find it, we ignore this match completely and let the search find it directly (rather than backMatching to find it)
		CorrelationWindow backMatchingWindow;

		try
		{
			if(centerOtherWindow == null) return null;

			float correlation = patchVolume.stretchCorrelations[iter];

			// set a penalty except for the last iterations
			float correlPenalty = 0.0f;
			if(iter < patchVolume.nIterations-1)
				correlPenalty = patchVolume.autoCorInc;

			// centerOtherWindow must be between bounding connections above and below and cannot cross
			// nextWindow is already between them, so move centerOtherWindow up or down to be between as well
			// new window selected must be of same type so move is +/- 4 index

			if(centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
			{
				if(StsPatchVolume.debugConnectCloseOnly) return null;
				//CorrelationWindow prevWindow = centerOtherWindow;
				while(centerOtherWindow != null && centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
				{
					int index = centerOtherWindow.windowIndex + 4;
					if(index >= otherTrace.nWindows)
						return null;
					//prevWindow = centerOtherWindow;
					centerOtherWindow = otherTrace.windows[index];
					if(centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
						return null;
				}
				//centerOtherWindow = newWindow.getClosestWindow(prevWindow, centerOtherWindow);
			}
			else if(centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
			{
				if(StsPatchVolume.debugConnectCloseOnly) return null;
				//CorrelationWindow prevWindow = centerOtherWindow;
				while(centerOtherWindow != null && centerOtherWindow.isBelowOrEqual(otherConnectWindowBelow))
				{
					int index = centerOtherWindow.windowIndex - 4;
					if(index < 0)
						return null;
					//prevWindow = centerOtherWindow;
					centerOtherWindow = otherTrace.windows[index];
					if(centerOtherWindow.isAboveOrEqual(otherConnectWindowAbove))
						return null;
				}
				//centerOtherWindow = newWindow.getClosestWindow(prevWindow, centerOtherWindow);
			}
			// centerOtherWindow is between above and below connection points on otherTrace, so compute correlation with window on this trace
			// this centerOtherWindow has already been determined to be the closest if their are two bracketing windows (@see
			if (newWindow.computeCorrelation(centerOtherWindow, CorrelationWindow.CENTER, correlPenalty) > correlation)
			{
				matchingWindow = centerOtherWindow;
				correlation = centerOtherWindow.stretchCorrelation;
			}

			if(StsPatchVolume.debugConnectCloseOnly)
				return matchingWindow;

			// try to make a match with otherWindow above centerOtherWindow
			centerOtherWindowIndex = centerOtherWindow.windowIndex;
			aboveWindowIndex = centerOtherWindowIndex - 4;
			if (aboveWindowIndex > otherConnectWindowAbove.windowIndex) // index must be below connectionAbove
			{
				aboveOtherWindow = otherTrace.windows[aboveWindowIndex];
				if (newWindow.computeCorrelation(aboveOtherWindow, CorrelationWindow.ABOVE, correlPenalty) > correlation)
				{
					matchingWindow = aboveOtherWindow;
					correlation = aboveOtherWindow.stretchCorrelation;
				}
			}

			belowWindowIndex = centerOtherWindowIndex + 4;
			if (belowWindowIndex < otherConnectWindowBelow.windowIndex)
			{
				belowOtherWindow = otherTrace.windows[belowWindowIndex];
				if (newWindow.computeCorrelation(belowOtherWindow, CorrelationWindow.BELOW, correlPenalty) > correlation)
				{
					matchingWindow = belowOtherWindow;
					// correlation = belowOtherWindow.stretchCorrelation;
				}
			}
			return matchingWindow;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(TracePoints.class, "connectWindows", e);
			return null;
		}
	}

	static CorrelationWindow getCenterOtherWindow(CorrelationWindow newWindow, TracePoints otherTrace, ConnectionList connectionList)
	{
		// center candidate window for matching
		CorrelationWindow centerOtherWindow;
		// bounding connections above and below which this connection cannot cross
		Connection connectionAbove, connectionBelow;
		// trace this new window is on
		TracePoints newTrace;
		// points on connections above and below
		PatchPoint newPointAbove, otherPointAbove, otherPointBelow;
		// offset from connectionAbove.point.traceIndex to the newWindow.centerPoint.traceIndex
		int newWindowIndexOffset;
		// index of candidate centerOtherWindow
		int centerOtherWindowIndex;
		// offset from the pointType at the parallel offset on the otherTrace from connectionAbove to the pointType we want on otherTrace
		int pointTypeOffset;
		// total offset including windowIndexOffset and pointTypeOffset
		int offset;

		connectionAbove = connectionList.connectionAbove;
		newTrace = newWindow.centerPoint.tracePoints;
		newPointAbove = connectionAbove.getWindowCenterPoint(newTrace);
		// newWindowIndexOffset is offset from connectionAbove to this newWindow.centerPoint
		newWindowIndexOffset = newWindow.windowIndex - newPointAbove.traceIndex;
		otherPointAbove = connectionAbove.getWindowCenterPoint(otherTrace);
		// compute offset to point on otherTrace which is parallel to connectionAbove and displaced down to newWindow
		// so the point is located newWindowIndexOffset below the otherPointAbove.traceIndex
		centerOtherWindowIndex = otherPointAbove.traceIndex + newWindowIndexOffset;
		if(centerOtherWindowIndex >= otherTrace.nWindows)
			centerOtherWindowIndex -= 4;
		centerOtherWindow = otherTrace.windows[centerOtherWindowIndex];
		// depending on the type here, we want to further offset from the type at this point to the desired type
		// which is equal to the numerical difference in type byte values
		// offset is always positive and is 0,1,2, or 3; if 0 or 1, use this offset; otherwise for 2, or 3, offset -2 or -1 respectively (offset-4)
		pointTypeOffset = getPointTypeDif(newWindow, centerOtherWindow);
		if (pointTypeOffset >= 2) pointTypeOffset -= 4;
		// the total offset is the shift from the connetionAbove (newWindowIndexOffset) plus the pointTypeOffset
		offset = newWindowIndexOffset + pointTypeOffset;
		if (offset < 1) offset += 4;
		centerOtherWindowIndex = otherPointAbove.traceIndex + offset;
		connectionBelow = connectionList.connectionBelow;
		otherPointBelow = connectionBelow.getWindowCenterPoint(otherTrace);

		if(centerOtherWindowIndex < otherPointBelow.traceIndex)
			return otherTrace.windows[centerOtherWindowIndex];

		while(centerOtherWindowIndex > otherPointAbove.traceIndex)
		{
			centerOtherWindowIndex -= 4;
			if(centerOtherWindowIndex < otherPointBelow.traceIndex)
				return otherTrace.windows[centerOtherWindowIndex];
		}
		return null;
	}

	static int getPointTypeDif(CorrelationWindow newWindow, CorrelationWindow otherWindow)
	{
		int typeDif = newWindow.centerPoint.pointType - otherWindow.centerPoint.pointType;
		if(typeDif < 0) typeDif += 4;
		return typeDif;
	}
	/**
	 * Given a newPatchPoint at newRow-newCol, which correlates with a prevPatchPoint at prevRow-prevCol which is possibly part of a patchGrid in the prevPatchGridsSet,
	 * combine these two points in the same patch.  The prevPatchPoint may be on the previous col (same row), or previous row (same col).
	 * If the previousPatchPoint is not part of an existing patchGrid (prevID == -1), then we will create a new patchGrid and add both points to it.
	 * If the previousPatchPoint is part of a patchGrid we will add the newPatchPoint to this patchGrid, unless the newPatchPoint already belongs to another patchGrid
	 * (this occurs when we first correlate with the previous column and find one patchGrid and also correlate with the previous row and find a different patchGrid).
	 */
	public Connection addPatchConnection(Connection connection, ConnectionList connectionList, boolean isRow)
	{
		StsPatchGrid patchGrid = null;
		CorrelationWindow window = connection.window;
		CorrelationWindow otherWindow = connection.otherWindow;
		float correlation = otherWindow.stretchCorrelation;
		// if (correlation < patchVolume.minLinkCorrel) return null;

		if(connectionList.connectionsCross(connection)) return null;

		PatchPoint newPatchPoint = window.centerPoint;
		PatchPoint otherPatchPoint = otherWindow.centerPoint;
		double distance = Math.abs(otherPatchPoint.slice - newPatchPoint.slice);

		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

		if(StsPatchVolume.debug && StsPatchGrid.debugPoint && (StsPatchGrid.doDebugPoint(newPatchPoint) || StsPatchGrid.doDebugPoint(otherPatchPoint)))
			StsException.systemDebug(this, "addPatchConnection", StsPatchVolume.iterLabel + " window " +
					newPatchPoint.toString() + " to " + otherPatchPoint.toString());

		// normally we can insert a new connectedPoint in the trace patchPointsList and split the connected interval;
		// but if we have cloned this new window and it is already connected, don't add/split the trace again
		//boolean splitIntervalOK = true;
		if (newPatchGrid == null)
		{
			if (otherPatchGrid == null) // prevPatchGrid doesn't exist, so create it and add otherPoint to it
			{
				patchGrid = StsPatchGrid.construct(patchVolume, newPatchPoint.getPointType(patchVolume.useFalseTypes));
				patchGrid.addPatchPoint(otherPatchPoint);
			}
			else // otherPatchGrid does exist, so use it
			{
				// if this newPatchPoint overlaps the otherPatchGrid, we can't add it;
				// So create a new patch and add a clone of the otherPatchPoint
				// try skipping for now...
				if (StsPatchVolume.debugCloneOK && otherPatchGrid.patchPointOverlaps(newPatchPoint)) // return null;
				{
					patchGrid = StsPatchGrid.construct(patchVolume, newPatchPoint.pointType);
					otherPatchPoint = otherPatchPoint.cloneAndClear();
					patchGrid.addPatchPoint(otherPatchPoint);
					//splitIntervalOK = false;
				}
				else // no overlap, so we will only need to add the newPatchPoint to it (below else)
					patchGrid = otherPatchGrid;
			}
			patchGrid.addPatchPoint(newPatchPoint);
		}
		else // newPatchGrid != null which means this window was just added to a patch from prevColTrace and the patchGrid would have been added to the rowGrids array
		{
			if (otherPatchGrid == null) // the otherPoint is not assigned to a patch; assign it to this one; don't add window to rowGrids unless it overlaps and addedGrid created
			{
				patchGrid = newPatchGrid;
				// otherPatchPoint doesn't have a patchGrid, but newPatchPoint does; try to add otherPatchPoint to newPatchGrid,
				// but if it overlaps, created an new patchGrid containing otherPatchPoint and a clone of newPatchPoint
				if (patchGrid.patchPointOverlaps(otherPatchPoint))
				{
					patchGrid = StsPatchGrid.construct(patchVolume, patchGrid.patchType);
					patchGrid.addPatchPoint(newPatchPoint.cloneAndClear());
				}
				patchGrid.addPatchPoint(otherPatchPoint);

				// patchGrid = patchGrid.checkAddPatchPoint(otherPatchPoint, newPatchPoint);
				//patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
				//checkAddPatchGridToRowGrids(patchGrid);
			}
			else if (otherPatchGrid.id == newPatchGrid.id) // otherPoint is already assigned to the same patch: addCorrelation
			{
				patchGrid = newPatchGrid;
				if (patchGrid == null) return null; // patchGrid not found; systemDebugError was printed
				//checkAddPatchGridToRowGrids(patchGrid);
				//patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correl);
			}
			// prevPoint and this window belong to different patches: merge newPatchGrid into prevPatchGrid and add connection
			// if we can't merge OK, then we create a new patch with newPoint and clone of connected otherPoint
			// cloned window is orphaned and won't be checked for additional connections; this will be done by the connected otherPoint.
			else
			{
				if (StsPatchGrid.mergePatchPointsOK(otherPatchGrid, newPatchGrid))
				{
					patchGrid = patchVolume.mergePatchGrids(otherPatchPoint, newPatchPoint);
					if (patchGrid == null)
						return null; // error occurred: systemError written in mergePatchGrids routine
				}
				// we can't merge grids, so add a cloned point of the newPatchPoint to the otherPatch
				// unless the otherPatch already has a window there; in this case, create a clone of the otherPatchPoint
				// and add it to a newPatchGrid with the newPoint
				// to anyone else via the trace search which doesn't see cloned points
				else if (StsPatchVolume.debugCloneOK)
				{
					if (!newPatchGrid.patchPointOverlaps(otherPatchPoint))
					{
						otherPatchPoint = otherPatchPoint.cloneAndClear();
						newPatchGrid.addPatchPoint(otherPatchPoint);
						patchGrid = newPatchGrid;
					}
					else if (!otherPatchGrid.patchPointOverlaps(newPatchPoint))
					{
						newPatchPoint = newPatchPoint.cloneAndClear();
						otherPatchGrid.addPatchPoint(newPatchPoint);
						patchGrid = otherPatchGrid;
					}
					else // each window overlaps the other grid, so we create a newGrid with clones of both
					{
						patchGrid = StsPatchGrid.construct(patchVolume, otherPatchGrid.patchType);
						newPatchPoint = newPatchPoint.cloneAndClear();
						patchGrid.addPatchPoint(newPatchPoint);
						otherPatchPoint = otherPatchPoint.cloneAndClear();
						patchGrid.addPatchPoint(otherPatchPoint);
						//splitIntervalOK = false;
					}
				}
			}
		}
		if (patchGrid == null) return null;

		checkResetClonedPoints(connection, otherPatchPoint, newPatchPoint, isRow);
		patchGrid.addCorrelation(connection, isRow);
		patchVolume.checkAddPatchGridToRowGrids(patchGrid);
		//if (splitIntervalOK)
		return connection;
		//else
		//	return null;
	}

	private void checkResetClonedPoints(Connection connection, PatchPoint otherPatchPoint, PatchPoint newPatchPoint, boolean isRow)
	{
		if (connection.otherWindow.centerPoint != otherPatchPoint) // otherPatchPoint may have been cloned
		{
			connection.otherWindow.centerPoint = otherPatchPoint;
			otherPatchPoint.setConnection(connection, isRow);
		}
		if (connection.window.centerPoint != newPatchPoint) // otherPatchPoint may have been cloned
		{
			connection.window.centerPoint = newPatchPoint;
			newPatchPoint.setConnection(connection, isRow);
		}
	}

	private CorrelationWindow findOtherNearestWindow(CorrelationWindow window, ConnectionList connectionList)
	{
		// first guess on window closest is at the same window array index
		// Use this as a start and search up and down for nearest
		int slice = window.centerPoint.slice;
		int otherIndex = getBoundedIndex(window.windowIndex);
		int otherIndexAbove = connectionList.getOtherPointIndexAbove();
		int otherIndexBelow = connectionList.getOtherPointIndexBelow();
		otherIndex = StsMath.limitBetweenExclusive(otherIndex, otherIndexAbove, otherIndexBelow);
		CorrelationWindow otherWindow = windows[otherIndex];

		int otherSlice = otherWindow.centerPoint.slice;
		byte pointType = StsTraceUtilities.coercedPointTypes[window.centerPoint.pointType];

		// find first nearest by searching forward and back
		CorrelationWindow windowBelow, windowAbove;
		Iterator<CorrelationWindow> forwardIterator = new WindowPointTypeForwardIterator(pointType, otherWindow, connectionList);
		windowBelow = forwardIterator.next();
		Iterator<CorrelationWindow> backwardIterator = new WindowPointTypeForwardIterator(pointType, otherWindow, connectionList);
		windowAbove = forwardIterator.next();

		if(windowAbove == null || windowAbove.centerPoint.slice < otherSlice)
		{
			if(windowBelow.centerPoint.slice <= otherSlice)
				return windowBelow;
			else
			   windowAbove = windowBelow;
		}
		else if(windowBelow == null)
		{
			if(windowAbove.centerPoint.slice >= otherSlice)
				return windowAbove;
			else
				windowAbove = windowBelow;
		}
		if (otherSlice < slice)
		{
			forwardIterator = new WindowPointTypeForwardIterator(pointType, otherWindow, connectionList);
			if(!forwardIterator.hasNext()) return null;
			CorrelationWindow nextWindow = forwardIterator.next();
			CorrelationWindow prevWindow;
			int nextDistance = nextWindow.centerPoint.slice - slice;
			while(forwardIterator.hasNext())
			{
				prevWindow = nextWindow;
				int prevDistance = nextDistance;
				nextWindow = forwardIterator.next();

				nextDistance = nextWindow.centerPoint.slice - slice;
				if(prevDistance <= 0 && nextDistance > 0)
				{
					if(-prevDistance < nextDistance)
						return prevWindow;
					else
						return nextWindow;
				}
			}
			return nextWindow;
		}
		else // otherSlice > slice
		{
			Iterator<CorrelationWindow> backwardsIterator = new WindowPointTypeBackwardIterator(pointType, otherWindow, connectionList);
			CorrelationWindow nextWindow = backwardsIterator.next();
			if(nextWindow == null) return null;
			CorrelationWindow prevWindow;
			int nextDistance = nextWindow.centerPoint.slice - slice;
			while(backwardsIterator.hasNext())
			{
				prevWindow = nextWindow;
				int prevDistance = nextDistance;
				nextWindow = backwardsIterator.next();

				nextDistance = nextWindow.centerPoint.slice - slice;
				if(prevDistance >= 0 && nextDistance < 0)
				{
					if(prevDistance > nextDistance)
						return prevWindow;
					else
						return nextWindow;
				}
			}
			return nextWindow;
		}
	}

	CorrelationWindow getWindowOfTypeBelow(CorrelationWindow window, ConnectionList connectionList)
	{
		Iterator<CorrelationWindow> forwardIterator = new WindowPointTypeForwardIterator(window.centerPoint.pointType, window, connectionList);
		return forwardIterator.next();
	}

	class WindowPointTypeForwardIterator implements Iterator<CorrelationWindow>
	{
		byte pointType;
		CorrelationWindow window;
		int otherIndexBelow;

		WindowPointTypeForwardIterator(byte pointType, CorrelationWindow window, ConnectionList connectionList)
		{
			this.pointType = pointType;
			this.window = window;
			otherIndexBelow = connectionList.getOtherPointIndexBelow();
			initialize(connectionList);
		}

		void initialize(ConnectionList connectionList)
		{
			if(window.centerPoint.pointType == pointType) return;
			window = getWindowBelowOfType(window, pointType, otherIndexBelow);
		}

		public boolean hasNext()
		{
			return window != null;
		}

		public CorrelationWindow next()
		{
			CorrelationWindow currentWindow = window;
			window = getWindowBelowOfType(window, pointType, otherIndexBelow);
			return currentWindow;
		}

		public void remove() {}
	}

	CorrelationWindow getWindowOfTypeAbove(CorrelationWindow window, ConnectionList connectionList)
	{
		Iterator<CorrelationWindow> backwardIterator = new WindowPointTypeBackwardIterator(window.centerPoint.pointType, window, connectionList);
		return backwardIterator.next();
	}

	CorrelationWindow getNextWindow(CorrelationWindow prevWindow)
	{
		int index = prevWindow.windowIndex + 1;
		if(index >= nWindows) return null;
		return windows[index];
	}

	class WindowPointTypeBackwardIterator implements Iterator<CorrelationWindow>
	{
		byte pointType;
		CorrelationWindow window;
		int otherIndexBelow;

		WindowPointTypeBackwardIterator(byte pointType, CorrelationWindow window, ConnectionList connectionList)
		{
			this.pointType = pointType;
			this.window = window;
			otherIndexBelow = connectionList.getOtherPointIndexBelow();
			if(window.centerPoint.pointType == pointType) return;
			this.window = getWindowAboveOfType(window, pointType, otherIndexBelow);
		}

		public boolean hasNext()
		{
			return window != null;
		}

		public CorrelationWindow next()
		{
			CorrelationWindow currentWindow = window;
			window = getWindowAboveOfType(window, pointType, otherIndexBelow);
			return currentWindow;
		}

		public void remove() {}
	}

	CorrelationWindow getWindowBelowOfType(CorrelationWindow window, byte pointType, int otherIndexBelow)
	{
		if(window == null) return null;
		int otherIndex = window.windowIndex;
		for(int i = otherIndex+1; i <= otherIndexBelow-1; i++)
			if(windows[i].centerPoint.pointType == pointType) return windows[i];
		return null;
	}

	CorrelationWindow getWindowAboveOfType(CorrelationWindow window, byte pointType, int otherIndexAbove)
	{
		int otherIndex = window.windowIndex;
		for(int i = otherIndex-1; i >= otherIndexAbove+1; i--)
			if(windows[i].centerPoint.pointType == pointType) return windows[i];
		return null;
	}

	int getBoundedIndex(int windowIndex)
	{
		return StsMath.minMax(windowIndex, 0, nWindows);
	}

	/** For this new window, we may have a new row and/or col connection or no connection.
	 *  If we have any new connection, then split our bounded connection interval at the new window
	 *  and move the interval down to this new interval. If the window was cloned, we need to use the
	 *  original window (window.clonedPoint) for these operations as it has the trace links for the
	 *  split and move operations.  If no connections, just move the interval down.
	 * @param window
	 * @param colConnection connection from the prevColPoint (same col, prev row) to this new window or its clone
	 * @param rowConnection connection from the prevRowPoint (same row, prev col) to this new window or its clone
	 */
	private void processNewConnections(CorrelationWindow window, Connection colConnection, Connection rowConnection)
	{
		try
		{
			// if we have a new row and/or col connection, split the connection interval
			// by inserting this connection into it
			if (colConnection != null)
			{
				colConnection = addPatchConnection(colConnection, colConnections, false);
				if(colConnection != null) colConnections.insert(colConnection);

			}
			if (rowConnection != null)
			{
				rowConnection = addPatchConnection(rowConnection, rowConnections, true);
				if(rowConnection != null) rowConnections.insert(rowConnection);
			}
			PatchPoint windowCenterPoint = window.centerPoint;
			if(windowCenterPoint.colConnection != null) colConnections.movePatchInterval(windowCenterPoint.colConnection);
			if(windowCenterPoint.rowConnection != null) rowConnections.movePatchInterval(windowCenterPoint.rowConnection);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "processNewConnections", e);
		}
	}
	/*
		private Connection checkAddRowConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
		{
			if (otherTrace == null) return null;
			if (window.centerPoint.rowConnection != null) return null;
			return window.checkAddConnection(otherTrace, minStretchCorrelation, true);
		}

		private Connection checkAddColConnection(CorrelationWindow window, TracePoints otherTrace, float minStretchCorrelation)
		{
			if (otherTrace == null) return null;
			if (window.centerPoint.colConnection != null) return null;
			return window.checkAddConnection(otherTrace, minStretchCorrelation, false);
		}

		Connection addConnection(boolean isRow, CorrelationWindow otherWindow, CorrelationWindow window, float correlation)
		{
			return new Connection(otherWindow, window, correlation);

		}
    */
	private void constructTraceWindows()
	{
		windows = new CorrelationWindow[nTracePatchPoints];
		nWindows = 0;
		PatchPoint prevPoint;
		PatchPoint point = null;
		PatchPoint nextPoint = tracePatchPoints[0];
		CorrelationWindow window;

		try
		{
			for (int n = 1; n < nTracePatchPoints - 1; n++)
			{
				prevPoint = point;
				point = nextPoint;
				nextPoint = tracePatchPoints[n];

				window = checkCreateWindow(prevPoint, point, nextPoint, nWindows);
				if (window == null) continue;
				if(nWindows != window.windowIndex)
					StsException.systemDebug(this, "constructTraceWindows", "Index out of sequence.");
				windows[nWindows] = window;
				nWindows++;
			}
			if (nWindows < nTracePatchPoints)
				windows = (CorrelationWindow[]) StsMath.trimArray(windows, nWindows);
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "constructTraceWindows", e);
		}
	}

	private CorrelationWindow checkCreateWindow(PatchPoint prevPoint, PatchPoint point, PatchPoint nextPoint, int windowIndex)
	{
		if(!arePointTypesOK(prevPoint, point, nextPoint)) return null;
		return new CorrelationWindow(prevPoint, point, nextPoint, windowIndex);
	}

	private boolean arePointTypesOK(PatchPoint prevPoint, PatchPoint point, PatchPoint nextPoint)
	{
		if(prevPoint != null)
		{
			if(nextPoint != null)
				return StsTraceUtilities.arePointTypesOK(prevPoint.pointType, point.pointType, nextPoint.pointType);
			else
				return StsTraceUtilities.arePointTypesAboveOK(prevPoint.pointType, point.pointType);
		}
		else // prevPoint == null
		{
			if(nextPoint != null)
				return StsTraceUtilities.arePointTypesBelowOK(point.pointType, nextPoint.pointType);
			else
				return false;
		}
	}

	private void reinitializeTraceIndices(TracePoints prevRowTrace, TracePoints prevColTrace)
	{
		if (this != null) reinitializeTraceIndexing();
		if (prevRowTrace != null) prevRowTrace.reinitializeTraceIndexing();
		if (prevColTrace != null) prevColTrace.reinitializeTraceIndexing();
	}

	void reinitializeTraceIndexing()
	{
		if(rowConnections != null) rowConnections.reinitializeTraceIndexing();
		if(colConnections != null) colConnections.reinitializeTraceIndexing();
	}

	/**
	 * currentPoint is the last currentPoint on this trace in the previous search operation, so is a good starting window for this search
	 * @param slice slice for which we want to find the nearest tracePoint
	 * @return the nearestTracePoint
	 */
	/*
		private PatchPoint nearestPatchPoint(int slice)
		{
			int distance;
			// if currentNearestPoint is above slice, search down for nearest
			if (currentPoint.slice < slice)
			{
				distance = slice - currentPoint.slice;
				PatchPoint point = currentPoint;

				for (int index = currentPoint.traceIndex + 1; index < nTracePatchPoints; index++)
				{
					PatchPoint lastPoint = point;
					int lastDistance = distance;
					point = tracePatchPoints[index];
					// if window is now below slice, then we have bracketed window: set and return currentPoint
					if (point.slice >= slice)
					{
						distance = point.slice - slice;
						if (distance < lastDistance)
							currentPoint = point;
						else
							currentPoint = lastPoint;
						return currentPoint;
					}
					else
						distance = slice - point.slice;
				}
				// didn't bracket, so slice is still below last window; return last window
				currentPoint = point;
			}
			// if currentNearestPoint is below slice, search up for nearest
			else if (currentPoint.slice > slice)
			{
				distance = currentPoint.slice - slice;
				PatchPoint point = currentPoint;

				for (int index = currentPoint.traceIndex - 1; index >= 0; index--)
				{
					PatchPoint lastPoint = point;
					int lastDistance = distance;
					point = tracePatchPoints[index];

					// if window is now above slice, then we have bracketed window: set and return currentPoint
					if (point.slice <= slice)
					{
						distance = slice - point.slice;
						if (distance < lastDistance)
							currentPoint = point;
						else
							currentPoint = lastPoint;
						return currentPoint;
					}
				}
				currentPoint = point;
			}
			return currentPoint;
		}

		private PatchPoint getOtherConnectedPatchPointAbove(boolean isRow)
		{
			return patchPointsList.getConnectionAbove(isRow).otherPoint;
		}

		private PatchPoint getOtherConnectedPatchPointBelow(boolean isRow)
		{
			return patchPointsList.getConnectionBelow(isRow).otherPoint;
		}
	*/
}

class PatchPoint implements Comparable<PatchPoint>, Cloneable
{
	/** trace this point is on */
	TracePoints tracePoints;
	float value;
	float z = StsParameters.nullValue;
	byte pointType;
	StsPatchGrid patchGrid;
	int slice;

	PatchPoint next, prev;
	/** connection from this tracePoint to the tracePoint on the adjacent trace at this row, col-1 (same row) */
	Connection rowConnection;
	/** connection from this tracePoint to the tracePoint on the adjacent trace at this row-1, col (same col) */
	Connection colConnection;
	/** index of this window in the trace containing it */
	int traceIndex;
	/** correl factor between this window and next in row. Note that rowConnection is from this window back. */
	float rowCorrel;
	/** correl factor between this window and next in col.  Note that colConnection is from this window back. */
	float colCorrel;
	/** cloned window for debugging.  Point this window was cloned from if cloned. */
	PatchPoint clonedPoint;

	/** first window above which has a connected patch */
//		PatchPoint connectionAbove = null;

	/** first window below which has a connected patch */
//		PatchPoint connectionBelow = null;

	PatchPoint()
	{
	}

	/** constructor for first and last links in doubly-linked list of PatchPoints */
	// PatchPoint(int traceIndex) { this.traceIndex = traceIndex; }

	/** constructor for first and last links in doubly-linked ConnectionList */
	PatchPoint(TracePoints trace, int slice)
	{
		this.tracePoints = trace;
		this.slice = slice;
		traceIndex = slice;
	}

	PatchPoint(TracePoints tracePoints, int row, int col, int slice, float z, float value, byte pointType, int traceIndex)
	{
		this.tracePoints = tracePoints;
		this.slice = slice;
		this.z = z;
		this.value = value;
		this.pointType = pointType;
		this.traceIndex = traceIndex;
	}

	PatchPoint(PatchPoint patchPoint, byte pointType, int traceIndex)
	{
		tracePoints = patchPoint.tracePoints;
		slice = patchPoint.slice;
		z = patchPoint.z;
		value = patchPoint.value;
		this.pointType = pointType;
		this.traceIndex = traceIndex;
	}


	PatchPoint(TracePoints tracePoints, int slice, int traceIndex)
	{
		this.tracePoints = tracePoints;
		this.slice = slice;
		this.traceIndex = traceIndex;
	}

	public boolean hasConnection()
	{
		return rowConnection != null || colConnection != null;
	}

	public boolean hasConnection(boolean isRow)
	{
		if (isRow)
			return rowConnection != null;
		else
			return colConnection != null;
	}

	public Connection getConnection(boolean isRow)
	{
		if (isRow) return rowConnection;
		else return colConnection;
	}

	public int compareTo(PatchPoint otherPoint)
	{
		if (slice > otherPoint.slice) return 1;
		else if (slice < otherPoint.slice) return -1;
		else return 0;
	}

	protected PatchPoint clone()
	{
		try
		{
			PatchPoint clonedPoint = (PatchPoint) super.clone();
			clonedPoint.clonedPoint = this;
			return clonedPoint;
		}
		catch (Exception e)
		{
			StsException.systemError(this, "clone");
			return null;
		}
	}

	/**
	 * A new point needs to be connected to a grid which already has a point at this location.
	 * So clone the window and clear any connection data.  This point will be added to the otherGrid
	 * defined by the otherPoint or to a new grid if there is no grid associated with the otherPoint.
	 * @return the cloned window
	 */
	protected PatchPoint cloneAndClear()
	{
		try
		{
			PatchPoint clonedPoint = (PatchPoint) super.clone();
			clonedPoint.clearConnectionData();
			clonedPoint.clonedPoint = this;
			return clonedPoint;
		}
		catch (Exception e)
		{
			StsException.systemError(this, "cloneAndClear");
			return null;
		}
	}

	void clearConnectionData()
	{
		rowConnection = null;
		colConnection = null;
		patchGrid = null;
		rowCorrel = 0.0f;
		colCorrel = 0.0f;
	}

	void setConnection(Connection connection, boolean isRow)
	{
		if(isRow)
			rowConnection = connection;
		else
			colConnection = connection;
	}

	Integer hashCode(int nVolumeCols)
	{
		return new Integer(getRow() * nVolumeCols + getCol());
	}

	int getSlice()
	{
		return slice;
	}

	int getID()
	{
		if (patchGrid == null) return -1;
		else return patchGrid.id;
	}

	int getIndex(int nVolumeCols)
	{
		return getCol() + getRow() * nVolumeCols;
	}

	public String toString()
	{
		if (clonedPoint != null)
			return pointToString() + " cloned from " + clonedPoint.patchToString();
		else
			return pointToString();
	}

	private String patchToString()
	{
		int id = -1;
		if (patchGrid != null) id = patchGrid.id;
		return "id " + id + " ";
	}

	private String pointToString()
	{
		return patchToString() + "r " + getRow() + " c " + getCol() + " s " + slice + " v " + value +
				" i " + traceIndex + " z " + z + " t " + StsTraceUtilities.typeStrings[pointType];
	}

	String nullOrToString(String string, PatchPoint patchPoint)
	{
		if (patchPoint == null) return " " + string + " null";
		else return " " + string + " " + patchPoint.toString();
	}

	StsPatchGrid getPatchGrid()
	{
		return patchGrid;
	}

	public void setPatchGrid(StsPatchGrid patchGrid)
	{
		this.patchGrid = patchGrid;
	}

	public byte getPointType(boolean useFalseTypes)
	{
		if (!useFalseTypes) return pointType;
		return StsTraceUtilities.coercedPointTypes[pointType];
	}

	public PatchPoint resetIndex(int index)
	{
		traceIndex = index;
		return this;
	}

	final protected int getRow()
	{
		return tracePoints.row;
	}

	final protected int getCol()
	{
		return tracePoints.col;
	}
}

/** A CorrelationWindow has a pointCenter and is bounded by trace points above and below of the appropriate point types.
 *  The window is essentially a half-wave.  Windows have the type of the center point.  A MAX window for example has
 *  a ZP (zero-plus) point above an a ZM (zero-minus) point below.  The skewness of the half-wave is considered in matching
 *  it to other windows so we retain for the window the minus and plus half-widths defined by the slice value difference.
 *  When matched with another windows, the stretchCorrelation is computed as the average of how much each side of the haf-wave
 *  has to be stretched to match the other.  A window with half-widths of -4 and +2 when matched with one with half-widths of
 *  -3 and +4 would have stretch ratios of .75 and 0.5 with an average of .667
 */
class CorrelationWindow implements Cloneable
{
	/** window in center of this window */
	PatchPoint centerPoint;
	/** windowType: BELOW if no window above center, ABOVE if no window below or CENTER */
	byte windowType;
	/** slice difference from center window to top window */
	int dSliceMinus;
	/** slice difference from bot window to center window */
	int dSlicePlus;
	/** correlation between this window and the connected window */
	float stretchCorrelation;
	/** index of this window in the windows array */
	int windowIndex;

	public static final byte CENTER = 0;
	public static final byte ABOVE = -1;
	public static final byte BELOW = 1;

	CorrelationWindow(PatchPoint pointAbove, PatchPoint centerPoint, PatchPoint pointBelow, int windowIndex)
	{
		this.centerPoint = centerPoint;
		this.windowIndex = windowIndex;
		if(pointAbove != null)
		{
			dSliceMinus = centerPoint.slice - pointAbove.slice;
			if (pointBelow != null)
			{
				dSlicePlus = pointBelow.slice - centerPoint.slice;
				windowType = CENTER;
			}
			else
			{
				dSlicePlus = dSliceMinus;
				windowType = ABOVE;
			}
		}
		else if (pointBelow != null) // pointAbove == null
		{
			dSlicePlus = pointBelow.slice - centerPoint.slice;
			dSliceMinus = dSlicePlus;
			windowType = BELOW;
		}
	}

	public CorrelationWindow clone()
	{
		try
		{
			CorrelationWindow window = (CorrelationWindow) super.clone();
			window.centerPoint = centerPoint.clone();
			return window;
		}
		catch (Exception e)
		{
			StsException.systemError(this, "clone");
			return null;
		}
	}

	boolean hasRowConnection() { return centerPoint.rowConnection != null; }

	boolean hasColConnection() { return centerPoint.colConnection != null; }

	boolean isAboveOrEqual(CorrelationWindow otherWindow)
	{
		return centerPoint.slice <= otherWindow.centerPoint.slice;
	}

	boolean isBelowOrEqual(CorrelationWindow otherWindow)
	{
		return centerPoint.slice >= otherWindow.centerPoint.slice;
	}

	/** check for closest of two windows where one must be above or equal to and the other must be below or equal to this window.
	 *  two windows can't be at same slice nor the order switched (above is below, below is above) as debug sanity checks
	 *
	 * @param windowAbove window above or equal in slice value to this window
	 * @param windowBelow window below or equal in slice value to this window
	 * @return
	 */
	CorrelationWindow getClosestWindow(CorrelationWindow windowAbove, CorrelationWindow windowBelow)
	{
		if(windowAbove == null)
			return windowBelow;
		else if(windowBelow == null)
			return windowAbove;
		int difAbove = centerPoint.slice - windowAbove.centerPoint.slice;
		int difBelow = windowBelow.centerPoint.slice - centerPoint.slice;
		if(StsPatchVolume.debug && (difAbove < 0 || difBelow < 0))
		{
			StsException.systemDebug(this, "getClosestWindow", "window not between windows above and below");
			return null;
		}
		if(StsPatchVolume.debug && (difAbove == 0 && difBelow == 0))
		{
//			StsException.systemDebug(this, "getClosestWindow", "other windows are the same, so can't be between");
			return null;
		}
		if(difAbove <= difBelow) return windowAbove;
		else					 return windowBelow;
	}

	/** check the various correlation measures and return the otherWindow if it matches; otherwise return null
	 *  store the correlation value in the otherWindow
	 * @param otherWindow correlation will be computed between this window and the otherWindow
	 * @param position
	 * @param correlPenalty
	 * @return correlation value (which is also stored in the otherWindow)
	 */
	float computeCorrelation(CorrelationWindow otherWindow, byte position, float correlPenalty)
	{
		// check correlation stretch
		float stretchCorrelation;
		if (windowType == BELOW)
			stretchCorrelation = computePlusStretchFactor(otherWindow);
		else if (windowType == ABOVE)
			stretchCorrelation = computeMinusStretchFactor(otherWindow);
		else
			stretchCorrelation = (computePlusStretchFactor(otherWindow) + computeMinusStretchFactor(otherWindow)) / 2;
		// store this correlation in both windows as well as returning it
		float totalPenalty = 0.0f;
		if(correlPenalty > 0.0f)
		{
			if(position != CENTER) totalPenalty = correlPenalty;
			if(StsTraceUtilities.isPointTypeFalse(centerPoint.pointType) || StsTraceUtilities.isPointTypeFalse(this.centerPoint.pointType))
				totalPenalty += correlPenalty;
		}
		stretchCorrelation -= totalPenalty;
		otherWindow.stretchCorrelation = stretchCorrelation;
		this.stretchCorrelation = stretchCorrelation;
		return stretchCorrelation;
	}
		/*
			private Connection checkAddConnection(TracePoints otherTrace, float minStretchCorrelation, boolean isRow)
			{
				CorrelationWindow otherMatchingWindow = findOtherMatchingWindow(otherTrace, isRow, minStretchCorrelation);
				if (otherMatchingWindow == null) return null;
				return addPatchConnection(otherMatchingWindow, isRow);
			}

			private ConnectionList getConnectionList(boolean isRow)
			{
				if(isRow) return rowConnections;
				else	  return colConnections;
			}

			boolean pointTypesMatch(CorrelationWindow otherWindow)
			{
				byte otherCenterType = otherWindow.centerPointType;
				if (centerPointType == otherCenterType) return true;
				if (!useFalseTypes) return false;

				byte centerType = StsTraceUtilities.coercedPointTypes[centerPointType];
				otherCenterType = StsTraceUtilities.coercedPointTypes[otherCenterType];
				return centerType == otherCenterType;
			}
        */
	/**
	 * check the above and below types to see that they match.
	 * We are assuming the centers have already been checked for matches
	 * @param window the window
	 * @param otherWindow otherWindow we are comparing it to
	 * @return true if centerTypes, and above and below types match
	 */
		/*
			boolean windowTypesMatch(CorrelationWindow window, CorrelationWindow otherWindow)
			{
				byte above = window.pointAbove.getPointType();
				byte below = window.pointBelow.getPointType();
				byte otherAbove = otherWindow.pointAbove.getPointType();
				byte otherBelow = otherWindow.pointBelow.getPointType();
				return above == otherAbove && below == otherBelow;
			}
        */
	/** if we have two windows with the exactly identical centerPoint, they must be equivalent if not equal windows. */
	boolean sameAs(CorrelationWindow otherWindow)
	{
		return otherWindow.centerPoint == centerPoint;
	}

	float computeMinusStretchFactor(CorrelationWindow otherWindow)
	{
		float minusStretchFactor = ((float) dSliceMinus) / otherWindow.dSliceMinus;
		if (minusStretchFactor > 1.0f)
			minusStretchFactor = 1 / minusStretchFactor;
		return minusStretchFactor;
	}

	float computePlusStretchFactor(CorrelationWindow otherWindow)
	{
		float plusStretchFactor = ((float) dSlicePlus) / otherWindow.dSlicePlus;
		if (plusStretchFactor > 1.0f)
			plusStretchFactor = 1 / plusStretchFactor;
		return plusStretchFactor;
	}
	/*
		boolean isCenterSliceOutsideWindow(int centerSlice)
		{
			return centerSlice < minSlice || centerSlice > maxSlice;
		}
	*/
	public String toString()
	{
		return " index " + windowIndex + " centerPoint: " + centerPoint.toString() + " correlation: " + stretchCorrelation;
	}
		/*
			private float computeStretchCorrelation(CorrelationWindow otherWindow)
			{
				if (otherWindow == null) return 0.0f;

				TracePoints traceOther = otherWindow.getTracePoints();
				int centerOther = otherWindow.centerSlice;
				int minOther = otherWindow.minSlice;
				int maxOther = otherWindow.maxSlice;

				// translate and stretch/shrink pointsOther z values so they line up with pointsNew

				int dSliceMinusOther = centerOther - minOther;
				int dSliceMinusNew = centerSlice - minSlice;
				float dSliceMinusOtherScalar = (float) dSliceMinusNew / dSliceMinusOther;
				// if(dzMinusOtherScalar < minStretchLimit || dzMinusOtherScalar > maxStretchLimit) return 0.0f;

				float minusStretchFactor = dSliceMinusOtherScalar;
				if (minusStretchFactor > 1.0f)
					minusStretchFactor = 1 / minusStretchFactor;

				int dSlicePlusOther = maxOther - centerOther;
				int dSlicePlusNew = maxSlice - centerSlice;
				float dSlicePlusOtherScalar = (float) dSlicePlusNew / dSlicePlusOther;
				// if(dzPlusOtherScalar < minStretchLimit || dzPlusOtherScalar > maxStretchLimit) return 0.0f;

				float plusStretchFactor = dSlicePlusOtherScalar;
				if (plusStretchFactor > 1.0f)
					plusStretchFactor = 1 / plusStretchFactor;

				return Math.min(minusStretchFactor, plusStretchFactor);
			}
			*/
}
/**
 * tracePoints has a series of row connections and col connections which are maintained during construction of this trace.
 * They will be removed after construction.  When checking on a new connection between this trace and the otherTrace
 * which is either row or col aligned with this trace, we bracket the search by connections above and below to prevent
 * crossing them.  On completion, if a new connection is created, it is added to to either trace.rowConnections or trace.colConnections.
 * These connection lists are currently double-linked, but could perhaps be only single-linked.
 */
class Connection
{
	Connection next, prev;
	/** connected window on this trace */
	CorrelationWindow window;
	/** connected window on other trace */
	CorrelationWindow otherWindow;
	/** average of slice values for two connected points; since connections can't cross or be identical, it will provide correct order */
	float sliceAvg;
	/** correlation between these two points; assigned to otherPoint location in either row or col direction */
	float correlation;

	Connection(CorrelationWindow otherWindow, CorrelationWindow newWindow)
	{
		this.window = newWindow;
		this.otherWindow = otherWindow;
		sliceAvg = (window.centerPoint.slice + otherWindow.centerPoint.slice)/2.0f;
		this.correlation = newWindow.stretchCorrelation;
	}

	PatchPoint getWindowCenterPoint(TracePoints tracePoints)
	{
		if(window.centerPoint.tracePoints == tracePoints) return window.centerPoint;
		else if(otherWindow.centerPoint.tracePoints == tracePoints) return otherWindow.centerPoint;
		StsException.systemDebug(this, "getWindowCenterPoint.  tracePoints ", tracePoints.toString() + " not found in " + toString());
		return null;
	}

	public String toString()
	{
		return " connection: " + window.toString() + " \n" + "     to other window " + otherWindow.toString() + " sliceAvg: " + sliceAvg + " correl: " + correlation;
	}
}
/**
 * doubly linked list of Connection[s].  There are two lists for each trace: rowConnections to the prev col on same row,
 * and colConnections to the prev row on the same col.  Connections in the list are in order, and must not cross or be identical.
 * Order is determined by the avg of the two slice values of the connected points (@see Connection).
 */
class ConnectionList
{
	/** first window in link list (connected to first actual window in list) */
	final Connection first;
	/** last window in link list (connected to last actual window in list) */
	final Connection last;
	/** last connected window in linked list just above current window */
	Connection connectionAbove;
	/** connected window just below connectionAbove in linked list */
	Connection connectionBelow;
	/** last connected window that was set; a convenient starting window for any search */
	// Connection currentConnection;

	/** Insert inactive row and col connections at the top and bottom of the connectionLists.
	 * New connections are added to these lists in order of the sliceAvg of the connection.
 	 * @param trace connection is from trace back to the otherTrace
	 * @param otherTrace trace connected which is either prevRow or prevCol
	 */
	ConnectionList(TracePoints trace, TracePoints otherTrace)
	{

		CorrelationWindow firstWindow = trace.windows[0].clone();
		firstWindow.windowIndex -= 1;
		firstWindow.centerPoint.slice -= 1;

		CorrelationWindow lastWindow = trace.windows[trace.nWindows-1].clone();
		lastWindow.windowIndex = trace.nWindows;
		lastWindow.centerPoint.slice += 1;

		CorrelationWindow firstOtherWindow = otherTrace.windows[0].clone();
		firstOtherWindow.windowIndex = -1;
		firstOtherWindow.centerPoint.slice -= 1;

		CorrelationWindow lastOtherWindow = otherTrace.windows[otherTrace.nWindows-1].clone();
		lastOtherWindow.windowIndex = otherTrace.nWindows;
		lastOtherWindow.centerPoint.slice += 1;

		first = new Connection(firstOtherWindow, firstWindow);
		last = new Connection(lastOtherWindow, lastWindow);
		first.next = last;
		last.prev = first;
		// currentConnection = first;
		connectionAbove = first;
		connectionBelow = last;
		connectionAbove.next = last;
		connectionBelow.prev = first;
	}

	void reinitializeTraceIndexing()
	{
		connectionAbove = first;
		connectionBelow = first.next;
		// currentConnection = first;
	}

	int getOtherPointIndexAbove()
	{
		return connectionAbove.otherWindow.windowIndex;
	}

	int getOtherPointIndexBelow()
	{
		return connectionBelow.otherWindow.windowIndex;
	}

	/** we have moved down to a new existing correlated patchPoint; set the interval to the one between this patchPoint and the window below */
	void movePatchInterval(Connection connectionAbove)
	{
		this.connectionAbove = connectionAbove;
		if(connectionAbove.next != null)
			connectionBelow = connectionAbove.next;
	}

	/**
	 * we are inserting this connectedPoint in an interval between connectionAbove and connectionBelow, either of which could be null
	 * meaning that it could be an openInterval with above and/or below undefined.  The interval (open or closed) is
	 * split into two subintervals and the current interval is set to the lower subinterval.
	 * @param connection between pointAbove and pointBelow where interval is to be split into two subIntervals.
	 */
	void insert(Connection connection)
	{
		if ( StsPatchVolume.debug && (connection == connectionAbove || connection == connectionBelow))
		{
			StsException.systemDebug(this, "insert", " connection " + connection.toString() + " same as " + connectionAbove.toString() + " or " + connectionBelow.toString());
			return;
		}
		connectionAbove.next = connection;
		connection.prev = connectionAbove;
		connection.next = connectionBelow;
		connectionBelow.prev = connection;

		// connectionAbove = connection;
		// connectionBelow = connection.next;
	}

	boolean connectionsCross(Connection connection)
	{
		if (!connectionsCross(connection, connectionAbove) && !connectionsCross(connection, connectionBelow))
			return false;

		if(StsPatchVolume.debug && StsPatchGrid.debugPoint && (StsPatchGrid.doDebugPoint(connection.window.centerPoint) || StsPatchGrid.doDebugPoint(connection.otherWindow.centerPoint)))
			StsException.systemDebug(this, "connectionCrosses", StsPatchVolume.iterLabel + connection.toString());

		return true;
	}

	boolean connectionsCross(Connection c1, Connection c2)
	{
		int crosses = StsMath.signProduct(c1.window.centerPoint.slice - c2.window.centerPoint.slice, c1.window.centerPoint.slice - c2.window.centerPoint.slice);
		return crosses < 0;
	}
}