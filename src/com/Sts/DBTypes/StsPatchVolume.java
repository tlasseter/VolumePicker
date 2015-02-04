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
	 * it is added to rowGrids.  If an existing grid is connected to a point in the row, it is added to rowGrid and removed from prevRowGrid.
	 * At the end of the row, grids still connected are in rowGrid, and disconnected ones are in prevRowGrids. These disconnected grids are
	 * added to gridList.  prevRowGrids is then set to rowGrids and rowGrids initialized for the next row.
	 */
	transient HashMap<Integer, StsPatchGrid> rowGrids = null;
	transient Iterator<StsPatchGrid> rowGridsIterator;

	/** prevRowGrids are the active patches in the previousRow; when making a connection to a point in the previous row, we look here for a patch. */
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
	/** pick point on adjoining trace cannot be more than this many wavelengths away from picked point */
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
	/** check that local loop doesn't cycle skip: if correls from prev row and prev col have same patch, allow only one connection point */
	transient boolean checkCycleSkips = false;
	/** For debugging: only run the cycle skip check if true */
	transient boolean cycleSkipOnly = false;
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

	public static final byte WINDOW_CENTERED = 0;
	public static final byte WINDOW_ABOVE = -1;
	public static final byte WINDOW_BELOW = 1;

	public final static byte POINT_ORIGINAL = StsTraceUtilities.POINT_ORIGINAL;

	public final static byte POINT_PLUS_ZERO_CROSSING = StsTraceUtilities.POINT_PLUS_ZERO_CROSSING;
	public final static byte POINT_MAXIMUM = StsTraceUtilities.POINT_MAXIMUM;
	public final static byte POINT_MINUS_ZERO_CROSSING = StsTraceUtilities.POINT_MINUS_ZERO_CROSSING;
	public final static byte POINT_MINIMUM = StsTraceUtilities.POINT_MINIMUM;

	public final static byte POINT_FALSE_MAXIMUM = StsTraceUtilities.POINT_FALSE_MAXIMUM;
	public final static byte POINT_FALSE_MINIMUM = StsTraceUtilities.POINT_FALSE_MINIMUM;

	public final static byte POINT_PLUS_FALSE_ZERO_CROSSING = StsTraceUtilities.POINT_PLUS_FALSE_ZERO_CROSSING;
	public final static byte POINT_MINUS_FALSE_ZERO_CROSSING = StsTraceUtilities.POINT_MINUS_FALSE_ZERO_CROSSING;

	public final static byte POINT_INTERPOLATED = StsTraceUtilities.POINT_INTERPOLATED;
	public final static byte POINT_FLAT_ZERO = StsTraceUtilities.POINT_FLAT_ZERO;

	public final static byte POINT_ANY = StsTraceUtilities.POINT_ANY;

	public final int largeInt = 99999999;

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
	static final boolean debug = false;
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
	static final int debugPatchID = StsPatchGrid.debugPatchID; // 71;
	static final boolean debugPatchGrid = debugPatchID != -1;
	static final boolean drawPatchBold = debugPatchGrid && StsPatchGrid.debugPatchGrid;
	static final boolean debugCloneOK = true;

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
		checkCycleSkips = pickPanel.checkCycleSkips;
		cycleSkipOnly = pickPanel.cycleSkipOnly;
		isIterative = pickPanel.isIterative;
		autoCorMax = pickPanel.autoCorMax;
		autoCorMin = pickPanel.autoCorMin;
		autoCorInc = pickPanel.autoCorInc;
		manualCorMin = pickPanel.manualCorMin;
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
		int croppedColMin = croppedBoundingBox.colMin;
		int croppedSliceMin = croppedBoundingBox.sliceMin;
		int nCroppedSlices = croppedBoundingBox.sliceMax - croppedSliceMin + 1;
		int nVolSlices = seismicVolume.nSlices;
		float[] traceValues = new float[nCroppedSlices];

		try
		{
			// row & col refer to the row and col in a croppedVolume over which picker is to run
			// volRow & volCol define the actual row and col in the volume (used only for reference)
			for (row = 0, volRow = croppedRowMin; row < nRows; row++, volRow++)
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
				for (col = 0, volCol = croppedColMin; col < nCols; col++, volCol++)
				{
					// StsException.systemDebug(this, "constructPatchVolume", "col loop, col: " + col);
					rowFloatBuffer.position(volCol * nSlices + croppedSliceMin);
					rowFloatBuffer.get(traceValues);

					TracePoints tracePoints = new TracePoints(row, col, traceValues, pickType, nSlices, croppedSliceMin, nVolSlices);
					rowTracePoints[col] = tracePoints;
					// prevColTracePoints are tracePoints in prev row & same col
					TracePoints prevColTracePoints = null;
					if (prevRowTracesPoints != null)
						prevColTracePoints = prevRowTracesPoints[col];

					// here we add the connected patchPoints
					tracePoints.addTracePatches(prevColTracePoints, prevRowTracePoints);

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

	class TracePoints
	{
		byte pickType;
		int row;
		int col;
		PatchPoint[] tracePatchPoints = new PatchPoint[0];
		int nTracePatchPoints;
		PatchPointLinkList patchPointsList;
		/** current tracePoint we are checking for connections */
		PatchPoint currentPoint;

		/** constructor for testing of link list (see main) */
		TracePoints(int nTracePatchPoints)
		{
			tracePatchPoints = new PatchPoint[nTracePatchPoints];
			for (int slice = 0; slice < nTracePatchPoints; slice++)
				tracePatchPoints[slice] = new PatchPoint(0, 0, slice, 0.0f, 0.0f, POINT_MAXIMUM, slice);
			initializePatchPointsList();
		}

		/**
		 * From the traceValues array, create an array of PatchPoint for each traceValue which qualifies as a pickType (all, min, max, +/- zero-crossing.
		 * This will be a sequential subset of the values array with nTracePatchPoints in this tracePatchPoints array
		 * @param row volume row of this trace
		 * @param col volume col of this trace
		 * @param traceValues original seismic values for this trace
		 * @param pickType pick type we are searching for (all, min, max, +/- zero-crossing; determines which events to pick
		 * @param nSlices number of total slices in this volume
		 * @param volSliceMin 0 or cropped volume slice min
		 * @param nVolSlices nSlices in volume or cropped volume
		 */
		TracePoints(int row, int col, float[] traceValues, byte pickType, int nSlices, int volSliceMin, int nVolSlices)
		{
			this.pickType = pickType;        // to match StsTraceUtilities.POINT...
			this.row = row;
			this.col = col;
			// tracePoints - uniform cubic interpolation of traceValues
			float[] tracePoints = StsTraceUtilities.computeCubicInterpolatedPoints(traceValues, nInterpolationIntervals);
			float z = zMin;
			if (tracePoints == null) return;
			int nTracePoints = tracePoints.length;
			tracePatchPoints = new PatchPoint[nTracePoints];
			volSliceMin *= nInterpolationIntervals;
			nVolSlices = (nVolSlices - 1) * nInterpolationIntervals + 1;
			int nTracePatchPoint = 0;
			for (int n = 0; n < nTracePoints; n++, z += interpolatedZInc)
			{
				byte tracePointType = StsTraceUtilities.getPointType(tracePoints, n);
				if (StsTraceUtilities.isMaxMinZeroOrFalseMaxMin(tracePointType))
					tracePatchPoints[nTracePatchPoint] = new PatchPoint(row, col, n, z, tracePoints[n], tracePointType, nTracePatchPoint++);
			}
			nTracePatchPoints = nTracePatchPoint;
			if (useFalseTypes)
				checkInsertFalseZeroCrossings(tracePoints);
			tracePatchPoints = (PatchPoint[]) StsMath.trimArray(tracePatchPoints, nTracePatchPoints);
			initializePatchPointsList();
		}

		/**
		 * creates correlated connections between this trace and traces at prev row & same col and prev col & same row
		 * @param prevColTrace prev trace in same col, prev row
		 * @param prevRowTrace prev trace in same row, prev col
		 */
		private void addTracePatches(TracePoints prevColTrace, TracePoints prevRowTrace)
		{
			if (prevRowTrace == null && prevColTrace == null) return;

			int nTracePatchPoints = tracePatchPoints.length;
			if (nTracePatchPoints == 0) return;

			// The patchPointsList for this trace is initialized with top and bot inactive connections to prev row and col traces.
			// These inactive connections are used to limit the search process.
			// a connection is from the first point in this patchPointsList to the first point in the corresponding row or col trace
			initializePatchPointsListConnections(prevColTrace, prevRowTrace);
			// If there are previous row and col traces, see if a point from each have the same patch;
			// if so, check whether a newPoint exists on this trace which correlates with these existing points.
			// If so, add it.
			if (checkCycleSkips && row > 0 && col > 0) addCyclePatches(prevColTrace, prevRowTrace, autoCorMin);
			// For each patchPoint, check previous traces for points which are same type and just above or below this new point
			// and find which of these two possible prev trace points has the best correlation.
			// If this correlation is above the minCorrelation, add this new point to the prev point patch.
			// If the new point already has a patch (because it correlated with one of the other of the 4 otherTraces, then the addPatchPointToPatch
			// will merge the two patches.

			if(cycleSkipOnly) return;

			for (int i = 0; i < nIterations; i++)
			{
				reinitializeTraceIndices(prevRowTrace, prevColTrace);
				float minStretchCorrelation = stretchCorrelations[i];
				for (int centerPointIndex = 0; centerPointIndex < nTracePatchPoints; centerPointIndex++)
				{
					PatchPoint centerPoint = tracePatchPoints[centerPointIndex];
					if (!StsTraceUtilities.isMaxMinOrZero(centerPoint.getPointType())) continue;
					CorrelationWindow window = constructCorrelationWindow(centerPoint);
					if (window != null)
						addTracePatch(window, prevColTrace, prevRowTrace, minStretchCorrelation);
				}
			}
		}

		/** If checkCycleSkips flag is set, this routine is the first pass for a trace and compares the three adjacent trace
		 *  in a loop with this trace at the upper right (e.g., if this trace is at 1,1 then traces at 0,0 & 0,1 & 1,0 form a loop).
		 *  These three other traces have already been processed, so have picked points and partial patches.
		 *  In our example 1,1 is the newTrace (this one), 0,1 is the prevColTrace and 1,0 is the prevRowTrace.  We arbitrarily
		 *  run down the prevColTrace looking for patches.  If a patch exists and also has a pint at the prevRowTrace, then we
		 *  want to check if a single point can be created at the newTrace.  We create temporary matchingWindows on the newTrace
		 *  correlating to windows on the prev row and col traces.  If these windows match and the correl coef is > minCorrelation,
		 *  we use these two new connections.  If they don't match, we have a potential cycle skip.  Taking the window which has the
		 *  highest correlation, compute the correlation coefficient between it and the other connection.  If this is ? minCorrelation,
		 *  we use them; otherwise we ignore and move on.
		 *
		 * @param prevColTrace trace on prev row, same col which forms one leg of square loop
		 * @param prevRowTrace trace on prev col, same row which forms the other leg of square loop
		 * @param minStretchCorrelation correlation required for acceptable connections
		 */
		private void addCyclePatches(TracePoints prevColTrace, TracePoints prevRowTrace, float minStretchCorrelation)
		{
			reinitializeTraceIndices(prevRowTrace, prevColTrace);
			PatchPointLinkList prevColConnectedPointList = prevColTrace.patchPointsList;
			PatchPoint prevColConnectedPoint;
			int prevRowTraceRow = prevRowTrace.row;
			int prevRowTraceCol = prevRowTrace.col;
			for (prevColConnectedPoint = prevColConnectedPointList.first; prevColConnectedPoint.next != null; prevColConnectedPoint = prevColConnectedPoint.next)
			{
				if (prevColConnectedPoint.patchGrid == null) continue;
				StsPatchGrid patchGrid = prevColConnectedPoint.patchGrid;
				PatchPoint prevRowConnectedPoint = patchGrid.getPatchPoint(prevRowTraceRow, prevRowTraceCol);
				if (prevRowConnectedPoint != null)
				{
					CorrelationWindow prevColWindow = prevColTrace.constructCorrelationWindow(prevColConnectedPoint);
					if (prevColWindow == null) continue;
					CorrelationWindow prevRowWindow = prevRowTrace.constructCorrelationWindow(prevRowConnectedPoint);
					if (prevRowWindow == null) continue;
					int newCenterSlice = (prevColConnectedPoint.slice + prevRowConnectedPoint.slice) / 2;
					PatchPoint newTraceCurrentPoint = nearestPatchPoint(newCenterSlice);
					CorrelationWindow matchingColWindow = prevColWindow.findNewMatchingWindow(this, false, minStretchCorrelation);
					if (matchingColWindow == null) continue;
					CorrelationWindow matchingRowWindow = prevRowWindow.findNewMatchingWindow(this, true, minStretchCorrelation);
					if (matchingRowWindow == null) continue;
					Connection colConnection, rowConnection;
					PatchPoint newPoint;
					// if prevRow and prevCol traces have points on same patch and they correlate to the same point on the new trace
					// add and process these connections
					if (matchingColWindow.sameAs(matchingRowWindow))
					{
						colConnection = matchingColWindow.addPatchConnection(prevColWindow, false);
						rowConnection = matchingRowWindow.addPatchConnection(prevRowWindow, true);
						newPoint = matchingColWindow.centerPoint;
						processNewConnections(newPoint, colConnection, rowConnection);
						continue;
					}
					// matching windows on this trace aren't the same for both prev row and col windows
					// so see if we can find one that best matches both
					// compare connecting prevColWindow with matchingRowWindow and prevRowWindow with matchingColWindow
					// take the best of these
					// colRowCorrelation is correlation from the prevCol to the matchingWindow for the prevRow
					// rowColCorrelation is correlation from the prevRow to thie matchingWindow for the prevCol
					float colRowCorrelation = matchingRowWindow.computeCorrelation(prevColWindow);
					float rowColCorrelation = matchingColWindow.computeCorrelation(prevRowWindow);
					// take the best of these and pair with corresponding other connection to create the two connections
					if (colRowCorrelation > rowColCorrelation && colRowCorrelation > minStretchCorrelation)
					{
						colConnection = matchingRowWindow.addPatchConnection(prevColWindow, false);
						rowConnection = matchingRowWindow.addPatchConnection(prevRowWindow, true);
						newPoint = matchingRowWindow.centerPoint;
						processNewConnections(newPoint, colConnection, rowConnection);
					}
					else if (rowColCorrelation > minStretchCorrelation)
					{
						rowConnection = matchingColWindow.addPatchConnection(prevRowWindow, true);
						colConnection = matchingColWindow.addPatchConnection(prevColWindow, false);
						newPoint = matchingColWindow.centerPoint;
						processNewConnections(newPoint, colConnection, rowConnection);
					}
				}
			}
		}

		private void addTracePatch(CorrelationWindow window, TracePoints prevColTrace, TracePoints prevRowTrace, float minStretchCorrelation)
		{
			PatchPoint newPoint = window.centerPoint;
			Connection rowConnection = null;
			Connection colConnection = null;

			if (prevColTrace != null)
				colConnection = checkAddColConnection(window, prevColTrace, minStretchCorrelation);
			if (prevRowTrace != null)
				rowConnection = checkAddRowConnection(window, prevRowTrace, minStretchCorrelation);

			processNewConnections(newPoint, colConnection, rowConnection);
		}

		protected CorrelationWindow constructCorrelationWindow(PatchPoint centerPoint)
		{
			try
			{
				return new CorrelationWindow(centerPoint);
			}
			catch (Exception e)
			{
				StsException.outputWarningException(this, "constructCorrelationWindow", e);
				return null;
			}
		}

		final boolean isMaximum(byte pointType)
		{
			return pointType == POINT_MAXIMUM || pointType == POINT_FALSE_MAXIMUM;
		}

		final boolean isMinimum(byte pointType)
		{
			return pointType == POINT_MINIMUM || pointType == POINT_FALSE_MINIMUM;
		}

		final boolean isZeroPlus(byte pointType)
		{
			return pointType == POINT_PLUS_ZERO_CROSSING || pointType == POINT_PLUS_FALSE_ZERO_CROSSING;
		}

		final boolean isZeroMinus(byte pointType)
		{
			return pointType == POINT_MINUS_ZERO_CROSSING || pointType == POINT_MINUS_FALSE_ZERO_CROSSING;
		}

		class CorrelationWindow
		{
			PatchPoint centerPoint;
			int centerPointIndex;
			PatchPoint pointAbove;
			PatchPoint pointBelow;
			/** slice value of center point */
			int centerSlice;
			/** slice value of top point */
			int minSlice;
			/** slice value of bot point */
			int maxSlice;
			/** slice difference from center point to top point */
			int dSliceMinus;
			/** slice difference from bot point to center point */
			int dSlicePlus;
			float stretchCorrelation;
			byte centerPointType;
			byte abovePointType;
			byte belowPointType;
			byte windowType;

			CorrelationWindow(PatchPoint centerPoint)
			{
				this.centerPoint = centerPoint;
				this.centerPointIndex = centerPoint.traceIndex;
				this.centerSlice = centerPoint.slice;
				centerPointType = centerPoint.getPointType();
				abovePointType = StsTraceUtilities.pointTypesBefore[centerPointType];
				belowPointType = StsTraceUtilities.pointTypesAfter[centerPointType];
				initialize();
			}

			boolean initialize()
			{
				try
				{
					int nTracePatchPoints = tracePatchPoints.length;
					if (nTracePatchPoints < 2) return false;

					if (centerPointIndex <= 0)
						return defineWindowBelow();
					else if (centerPointIndex == nTracePatchPoints - 1)
						return defineWindowAbove();
					else
						return defineWindowCentered();
				}
				catch (Exception e)
				{
					StsException.outputWarningException(this, "initialize", e);
					return false;
				}
			}

			private TracePoints getTracePoints()
			{
				return TracePoints.this;
			}

			boolean defineWindowBelow()
			{
				windowType = WINDOW_BELOW;
				pointAbove = centerPoint;
				pointBelow = getTracePatchPointBelow(centerPointIndex);
				if (pointBelow == null) return false;
				maxSlice = pointBelow.slice;
				dSlicePlus = maxSlice - centerSlice;
				dSliceMinus = dSlicePlus;
				minSlice = centerSlice - dSliceMinus;
				return true;
			}

			boolean defineWindowAbove()
			{
				windowType = WINDOW_ABOVE;
				pointBelow = centerPoint;
				pointAbove = getTracePatchPointAbove(centerPointIndex);
				if (pointAbove == null) return false;
				minSlice = pointAbove.slice;
				dSliceMinus = centerSlice - minSlice;
				dSlicePlus = dSliceMinus;
				maxSlice = centerSlice + dSlicePlus;
				return true;
			}

			boolean defineWindowCentered()
			{
				windowType = WINDOW_CENTERED;
				pointAbove = getTracePatchPointAbove(centerPointIndex);
				pointBelow = getTracePatchPointBelow(centerPointIndex);
				if (pointAbove == null && pointBelow == null) return false;
				if (pointAbove == null) return defineWindowBelow();
				if (pointBelow == null) return defineWindowAbove();
				minSlice = pointAbove.slice;
				maxSlice = pointBelow.slice;
				dSliceMinus = centerSlice - minSlice;
				dSlicePlus = maxSlice - centerSlice;
				return true;
			}

			PatchPoint getTracePatchPointBelow(int centerPointIndex)
			{
				for (int index = centerPointIndex + 1; index < nTracePatchPoints; index++)
					if (tracePatchPoints[index].getPointType() == belowPointType)
						return tracePatchPoints[index];
				return null;
			}

			PatchPoint getTracePatchPointAbove(int centerPointIndex)
			{
				for (int index = centerPointIndex - 1; index >= 0; index--)
					if (tracePatchPoints[index].getPointType() == abovePointType)
						return tracePatchPoints[index];
				return null;
			}

			/**
			 * Determine interval on otherTrace between two correlated patchPoints which bracket the centerSlice on thisTrace.
			 * Search up and down from centerSlice on otherTrace for anywhere from two to four possible matchingWindows.
			 * If a possible match is found within "near" distance from centerSlice, search for a second in that direction;
			 * otherwise just check the one.
			 * @param otherTrace the trace on which we want to find matching windows
			 * @param minCorrelation minimum value used in matching windows @see matches()
			 * @return best correlation window or null if none
			 */

			static final int near = 3;

			CorrelationWindow findOtherMatchingWindow(TracePoints otherTrace, boolean isRow, float minCorrelation)
			{
				// this is sliceIndex for centerPoint of this window for which we want to find a matchingWindow on otherTrace
				int centerSlice = centerPoint.slice;
				// sliceIndex for correlated point just above centerSlice on otherTrace otherTrace; search for matches down from here to correlatedPoint below
				PatchPoint otherCenterPoint = otherTrace.nearestPatchPoint(centerSlice);
				int otherCenterIndex = otherCenterPoint.traceIndex;

				PatchPoint otherConnectedPatchPointAbove = getOtherConnectedPatchPointAbove(isRow);
				PatchPoint otherConnectedPatchPointBelow = getOtherConnectedPatchPointBelow(isRow);
				int otherConnectedPointIndexAbove = otherConnectedPatchPointAbove.traceIndex;
				int otherConnectedPointIndexBelow = otherConnectedPatchPointBelow.traceIndex;
				otherCenterIndex = StsMath.limitBetweenExclusive(otherCenterIndex, otherConnectedPointIndexAbove, otherConnectedPointIndexBelow);
				if (otherCenterIndex == StsParameters.nullInteger) return null;

				return findMatchingWindow(otherTrace, otherCenterPoint, otherConnectedPointIndexAbove, otherConnectedPointIndexBelow, minCorrelation);
			}

			/**
			 * here we are looking for a match window on a new trace where thisTrace is the otherTrace.
			 * So we need to find the connectedPoints above and below the centerPoint on this newTrace which limit it.
			 * @param newTrace trace
			 * @param isRow
			 * @param minCorrelation
			 * @return
			 */
			CorrelationWindow findNewMatchingWindow(TracePoints newTrace, boolean isRow, float minCorrelation)
			{
				// this is sliceIndex for centerPoint of this otherWindow for which we want to find a matchingWindow on newTrace
				int otherCenterSlice = centerPoint.slice;
				// sliceIndex for correlated point just above centerSlice on trace; search for matches down from here to correlatedPoint below
				PatchPoint newStartPoint = newTrace.nearestPatchPoint(otherCenterSlice);

				PatchPoint newConnectedPatchPointAbove = newTrace.patchPointsList.getConnectedPatchPointAbove(newStartPoint, isRow);
				PatchPoint newConnectedPatchPointBelow = newConnectedPatchPointAbove.next;
				int newConnectedPointIndexAbove = newConnectedPatchPointAbove.traceIndex;
				int newConnectedPointIndexBelow = newConnectedPatchPointBelow.traceIndex;
				int newCenterIndex = newStartPoint.traceIndex;
				newCenterIndex = StsMath.limitBetweenExclusive(newCenterIndex, newConnectedPointIndexAbove, newConnectedPointIndexBelow);
				if (newCenterIndex == StsParameters.nullInteger) return null;
				PatchPoint newCenterPoint = newTrace.tracePatchPoints[newCenterIndex];
				return findMatchingWindow(newTrace, newCenterPoint, newConnectedPointIndexAbove, newConnectedPointIndexBelow, minCorrelation);
			}

			/**
			 * Find matchingWindow on otherTrace which matches this window.
			 * First determine interval on otherTrace between two correlated patchPoints which bracket the input startPoint on otherTrace.
			 * Search up and down from startPoint.traceIndex on otherTrace for anywhere from two to four possible matchingWindows.
			 * If a possible match is found within "near" distance from startSlice, search for a second in that direction;
			 * otherwise just check the one.
			 * @param otherTrace trace we want to search for the best matching window to this one
			 * @param otherStartPoint point on otherTrace which is starting point to search up and down for matching windows
			 * @param otherConnectedPointIndexAbove tracePatchPoints index which is upper limit of search (last connectedPoint above: we don't want connection to cross it)
			 * @param otherConnectedPointIndexBelow tracePatchPoints index which is lower limit of search (last connectedPoint below: we don't want connection to cross it)
			 * @param minCorrelation minimum value used in matching windows @see matches()
			 * @return best correlation window or null if none
			 */
			CorrelationWindow findMatchingWindow(TracePoints otherTrace, PatchPoint otherStartPoint, int otherConnectedPointIndexAbove, int otherConnectedPointIndexBelow, float minCorrelation)
			{
				// sliceIndex for correlated point just above centerSlice on otherTrace otherTrace; search for matches down from here to correlatedPoint below
				int startIndex = otherStartPoint.traceIndex;
				int startSlice = otherStartPoint.slice;

				CorrelationWindow matchingWindow = null; // this will be best matching window in otherTrace

				int offset = 0;
				int index;
				PatchPoint otherTraceCenterPoint;
				CorrelationWindow window;

				CorrelationWindow[] correlationWindows;
				int nCorrelationWindows = 0;

				if (debugCorrelationWindows)
					correlationWindows = new CorrelationWindow[4];

				try
				{
					// search down for one or two windows (two if first is within "near" slices of center
					for (index = startIndex; index < otherConnectedPointIndexBelow; index++, offset++)
					{
						otherTraceCenterPoint = otherTrace.tracePatchPoints[index];
						if (otherTraceCenterPoint.slice > maxSlice) break;
						if (!pointTypesMatch(centerPointType, otherTraceCenterPoint.getPointType())) continue;
						window = checkCreateOtherMatchingWindow(otherTrace, otherTraceCenterPoint, centerPointType);
						if (window == null) continue;

						if (debugCorrelationWindows)
							correlationWindows[nCorrelationWindows++] = window;
						// to be accepted, the window correlation must be the best of all and greater than the minCorrelation
						if (window.stretchCorrelation > minCorrelation)
						{
							matchingWindow = window;
							minCorrelation = window.stretchCorrelation;
						}

						if (!searchForMultipleWindowMatches || offset > near)
							break; // if we found one inside near, try for a second but once outside near, break give us one or two windows
					}
					offset = 0;
					// search up for one or two windows (two if first is within "near" slices of center
					for (index = startIndex - 1; index > otherConnectedPointIndexAbove; index--, offset++)
					{
						otherTraceCenterPoint = otherTrace.tracePatchPoints[index];
						if (otherTraceCenterPoint.slice < minSlice) break;
						if (centerPointType != otherTraceCenterPoint.getPointType()) continue;
						window = checkCreateOtherMatchingWindow(otherTrace, otherTraceCenterPoint, centerPointType);
						if (window == null) continue;

						if (debugCorrelationWindows)
							correlationWindows[nCorrelationWindows++] = window;
						if (window.stretchCorrelation > minCorrelation)
						{
							matchingWindow = window;
							minCorrelation = window.stretchCorrelation;
						}
						if (!searchForMultipleWindowMatches || offset > near)
							break;
					}
					if (debugCorrelationWindows && nCorrelationWindows > 1)
					{
						System.out.println("findNewMatchingWindow: " + nCorrelationWindows + " windows");
						for (int n = 0; n < nCorrelationWindows; n++)
							System.out.println(" window: " + n + correlationWindows[n].toString());
					}
					return matchingWindow;
				}
				catch (Exception e)
				{
					StsException.outputWarningException(this, "findNewMatchingWindow", e);
					return null;
				}
			}

			CorrelationWindow checkCreateOtherMatchingWindow(TracePoints otherTrace, PatchPoint otherTraceCenterPoint, byte centerPointType)
			{
				CorrelationWindow otherWindow = otherTrace.constructCorrelationWindow(otherTraceCenterPoint);
				if (otherWindow == null) return null;
				if (!windowTypesMatch(this, otherWindow)) return null;
				if (otherWindow.isCenterSliceOutsideWindow(centerSlice)) return null;
				computeCorrelation(otherWindow);
				return otherWindow;
			}

			/** check the various correlation measures and return the otherWindow if it matches; otherwise return null
			 *  store the correlation value in the otherWindow
			 * @param otherWindow correlation will be computed between this window and the otherWindow
			 * @return correlation value (which is also stored in the otherWindow)
			 */
			float computeCorrelation(CorrelationWindow otherWindow)
			{
				// check correlation stretch
				float stretchCorrelation;
				if (windowType == WINDOW_BELOW)
					stretchCorrelation = computePlusStretchFactor(otherWindow);
				else if (windowType == WINDOW_ABOVE)
					stretchCorrelation = computeMinusStretchFactor(otherWindow);
				else
					stretchCorrelation = (computePlusStretchFactor(otherWindow) + computeMinusStretchFactor(otherWindow)) / 2;
				// store this correlation in both windows as well as returning it
				otherWindow.stretchCorrelation = stretchCorrelation;
				this.stretchCorrelation = stretchCorrelation;
				return stretchCorrelation;
			}

			private Connection checkAddConnection(TracePoints otherTrace, float minStretchCorrelation, boolean isRow)
			{
				CorrelationWindow otherMatchingWindow = findOtherMatchingWindow(otherTrace, isRow, minStretchCorrelation);
				if (otherMatchingWindow == null) return null;
				return addPatchConnection(otherMatchingWindow, isRow);
			}

			/**
			 * Given a newPatchPoint at newRow-newCol, which correlates with a prevPatchPoint at prevRow-prevCol which is possibly part of a patchGrid in the prevPatchGridsSet,
			 * combine these two points in the same patch.  The prevPatchPoint may be on the previous col (same row), or previous row (same col).
			 * If the previousPatchPoint is not part of an existing patchGrid (prevID == -1), then we will create a new patchGrid and add both points to it.
			 * If the previousPatchPoint is part of a patchGrid we will add the newPatchPoint to this patchGrid, unless the newPatchPoint already belongs to another patchGrid
			 * (this occurs when we first correlate with the previous column and find one patchGrid and also correlate with the previous row and find a different patchGrid).
			 */
			public Connection addPatchConnection(CorrelationWindow otherWindow, boolean isRow)
			{
				StsPatchGrid patchGrid;
				float correlation = otherWindow.stretchCorrelation;
				if (correlation < minLinkCorrel) return null;

				PatchPoint newPatchPoint = centerPoint;
				PatchPoint otherPatchPoint = otherWindow.centerPoint;
				double distance = Math.abs(otherPatchPoint.slice - newPatchPoint.slice);

				StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
				StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();

				// normally we can insert a new connectedPoint in the trace patchPointsList and split the connected interval;
				// but if we have cloned this new point and it is already connected, don't add/split the trace again
				//boolean splitIntervalOK = true;
				if (newPatchGrid == null)
				{
					if (otherPatchGrid == null) // prevPatchGrid doesn't exist, so create it and add otherPoint to it
					{
						patchGrid = StsPatchGrid.construct(StsPatchVolume.this, newPatchPoint.getPointType());
						patchGrid.addPatchPoint(otherPatchPoint);
					}
					else // otherPatchGrid does exist, so use it
					{
						// if this newPatchPoint overlaps the otherPatchGrid, we can't add it;
						// So create a new patch and add a clone of the otherPatchPoint
						if (otherPatchGrid.patchPointOverlaps(newPatchPoint))
						{
							patchGrid = StsPatchGrid.construct(StsPatchVolume.this, newPatchPoint.pointType);
							otherPatchPoint = otherPatchPoint.clone();
							patchGrid.addPatchPoint(otherPatchPoint);
							//splitIntervalOK = false;
						}
						else // no overlap, so we will only need to add the newPatchPoint to it (below else)
							patchGrid = otherPatchGrid;
					}
					patchGrid.addPatchPoint(newPatchPoint);
				}
				else // id != -1 means this point was just added to a patch from prev connection and the patchGrid would have been added to the rowGrids array
				{
					if (otherPatchGrid == null) // the otherPoint is not assigned to a patch; assign it to this one; don't add point to rowGrids unless it overlaps and addedGrid created
					{
						patchGrid = newPatchGrid;
						if (patchGrid == null) return null; // patchGrid not found; systemDebugError was printed
						// otherPatchPoint doesn't have a patchGrid, but newPatchPoint does; try to add otherPatchPoint to newPatchGrid,
						// but if it overlaps, created an addedPatchGrid containing otherPatchPoint and a clone of newPatchPoint
						// return this addedPatchGrid or patchGrid (if no new grid added)
						patchGrid = patchGrid.checkAddPatchPoint(otherPatchPoint, newPatchPoint);
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
					// prevPoint and this point belong to different patches: merge newPatchGrid into prevPatchGrid and add connection
					// if we can't merge OK, then we create a new patch with newPoint and clone of connected otherPoint
					// cloned point is orphaned and won't be checked for additional connections; this will be done by the connected otherPoint.
					else
					{
						if (StsPatchGrid.mergePatchPointsOK(otherPatchGrid, newPatchGrid))
						{
							patchGrid = mergePatchGrids(otherPatchPoint, newPatchPoint);
							if (patchGrid == null)
								return null; // error occurred: systemError written in mergePatchGrids routine
						}
						else
						{
							// we cant merge grids, so add a cloned point of the newPatchPoint to the otherPatch
							// unless the otherPatch already has a point there; in this case, create a clone of the otherPatchPoint
							// and add it to a newPatchGrid with the newPoint
							// to anyone else via the trace search which doesn't see cloned points
							if (debugCloneOK)
							{
								if (!otherPatchGrid.patchPointOverlaps(newPatchPoint))
								{
									newPatchPoint = newPatchPoint.clone();
									otherPatchGrid.addPatchPoint(newPatchPoint);
									patchGrid = otherPatchGrid;
								}
								else if (!newPatchGrid.patchPointOverlaps(otherPatchPoint))
								{
									otherPatchPoint = otherPatchPoint.clone();
									newPatchGrid.addPatchPoint(otherPatchPoint);
									patchGrid = newPatchGrid;
								}
								else // each point overlaps the other grid, so we create a newGrid with clones of both
								{
									patchGrid = StsPatchGrid.construct(StsPatchVolume.this, otherPatchGrid.patchType);
									newPatchPoint = newPatchPoint.clone();
									patchGrid.addPatchPoint(newPatchPoint);
									otherPatchPoint = otherPatchPoint.clone();
									patchGrid.addPatchPoint(otherPatchPoint);
									//splitIntervalOK = false;
								}
							}
						}
					}
				}

				if (patchGrid == null) return null;

				patchGrid.addCorrelation(otherPatchPoint, newPatchPoint, correlation);
				Connection connection = addConnection(isRow, otherPatchPoint, newPatchPoint, correlation);
				checkAddPatchGridToRowGrids(patchGrid);
				//if (splitIntervalOK)
					return connection;
				//else
				//	return null;
			}

			boolean pointTypesMatch(byte centerType, byte otherCenterType)
			{
				if (centerType == otherCenterType) return true;

				if (!useFalseTypes) return false;

				centerType = StsTraceUtilities.coercedPointTypes[centerType];
				otherCenterType = StsTraceUtilities.coercedPointTypes[otherCenterType];
				return centerType == otherCenterType;
			}

			/**
			 * check the above and below types to see that they match.
			 * We are assuming the centers have already been checked for matches
			 * @param window the window
			 * @param otherWindow otherWindow we are comparing it to
			 * @return true if centerTypes, and above and below types match
			 */
			boolean windowTypesMatch(CorrelationWindow window, CorrelationWindow otherWindow)
			{
				byte above = window.pointAbove.getPointType();
				byte below = window.pointBelow.getPointType();
				byte otherAbove = otherWindow.pointAbove.getPointType();
				byte otherBelow = otherWindow.pointBelow.getPointType();
				return above == otherAbove && below == otherBelow;
			}

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

			boolean isCenterSliceOutsideWindow(int centerSlice)
			{
				return centerSlice < minSlice || centerSlice > maxSlice;
			}

			public String toString()
			{
				return "centerPoint: " + centerPoint.toString() + " correlation: " + stretchCorrelation;
			}

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
		}

		/** For this new point, we may have a new row and/or col connection or no connection.
		 *  If we have any new connection, then split our bounded connection interval at the new point
		 *  and move the interval down to this new interval. If the point was cloned, we need to use the
		 *  original point (point.clonedPoint) for these operations as it has the trace links for the
		 *  split and move operations.  If no connections, just move the interval down.
		 * @param newPoint new point which is at the end of the two row and/or col connections; connection may have a clone of it
		 * @param colConnection connection from the prevColPoint (same col, prev row) to this new point or its clone
		 * @param rowConnection connection from the prevRowPoint (same row, prev col) to this new point or its clone
		 */
		private void processNewConnections(PatchPoint newPoint, Connection colConnection, Connection rowConnection)
		{
			if (colConnection != null)
			{
				splitPatchInterval(colConnection.point);
				movePatchInterval(colConnection.point);
			}
			else if (rowConnection != null)
			{
				splitPatchInterval(rowConnection.point);
				movePatchInterval(rowConnection.point);
			}
			else
				movePatchInterval(newPoint);
		}

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

		Connection addConnection(boolean isRow, PatchPoint otherPatchPoint, PatchPoint patchPoint, float correlation)
		{
			Connection connection = new Connection(otherPatchPoint, patchPoint, correlation);
			if (isRow) patchPoint.rowConnection = connection;
			else patchPoint.colConnection = connection;
			return connection;
		}

		private void checkInsertFalseZeroCrossings(float[] tracePoints)
		{
			PatchPoint[] newTracePatchPoints = new PatchPoint[tracePatchPoints.length];

			PatchPoint point = tracePatchPoints[0];
			byte pointType = point.pointType;
			PatchPoint nextPoint = tracePatchPoints[1];
			byte nextPointType = nextPoint.pointType;

			int nn = 0; // index of new array with new points added
			newTracePatchPoints[nn++] = point;
			int row = point.row;
			int col = point.col;
			float z = point.z;
			for (int n = 1; n < nTracePatchPoints - 1; n++, z += interpolatedZInc)
			{
				PatchPoint prevPoint = point;
				byte prevPointType = pointType;
				point = nextPoint;
				pointType = nextPointType;
				nextPoint = tracePatchPoints[n + 1];
				nextPointType = nextPoint.pointType;
				if (isMaximum(pointType))
				{
					if (prevPointType == POINT_FALSE_MINIMUM)
					{
						float sliceDif = point.slice - prevPoint.slice;
						if (sliceDif > 2)
						{
							float floatNewSlice = sliceDif / 2 + prevPoint.slice;
							int newSlice = intervalRound(floatNewSlice, prevPoint.slice);
							float newZ = prevPoint.z + interpolatedZInc * (sliceDif / 2);
							newTracePatchPoints[nn] = new PatchPoint(row, col, newSlice, newZ, tracePoints[newSlice], POINT_PLUS_FALSE_ZERO_CROSSING, nn);
							nn++;
						}
					}

					newTracePatchPoints[nn] = point.resetIndex(nn);
					nn++;

					if (nextPointType == POINT_FALSE_MINIMUM)
					{
						float sliceDif = nextPoint.slice - point.slice;
						if (sliceDif > 2)
						{
							float floatNewSlice = sliceDif / 2 + point.slice;
							int newSlice = intervalRound(floatNewSlice, point.slice);
							float newZ = point.z + interpolatedZInc * (sliceDif / 2);
							newTracePatchPoints[nn] = new PatchPoint(row, col, newSlice, newZ, tracePoints[newSlice], POINT_MINUS_FALSE_ZERO_CROSSING, nn);
							nn++;
						}
					}
				}
				else if (isMinimum(pointType))
				{
					if (prevPointType == POINT_FALSE_MAXIMUM)
					{
						float sliceDif = point.slice - prevPoint.slice;
						if (sliceDif > 2)
						{
							float floatNewSlice = sliceDif / 2 + prevPoint.slice;
							int newSlice = intervalRound(floatNewSlice, prevPoint.slice);
							float newZ = prevPoint.z + interpolatedZInc * (sliceDif / 2);
							newTracePatchPoints[nn] = new PatchPoint(row, col, newSlice, newZ, tracePoints[newSlice], POINT_MINUS_FALSE_ZERO_CROSSING, nn);
							nn++;
						}
					}

					newTracePatchPoints[nn] = point.resetIndex(nn);
					nn++;

					if (nextPointType == POINT_FALSE_MAXIMUM)
					{
						float sliceDif = nextPoint.slice - point.slice;
						if (sliceDif > 2)
						{
							float floatNewSlice = sliceDif / 2 + point.slice;
							int newSlice = intervalRound(floatNewSlice, point.slice);
							float newZ = point.z + interpolatedZInc * (sliceDif / 2);
							newTracePatchPoints[nn] = new PatchPoint(row, col, newSlice, newZ, tracePoints[newSlice], POINT_PLUS_FALSE_ZERO_CROSSING, nn);
							nn++;
						}
					}
				}
				else
				{
					newTracePatchPoints[nn] = point.resetIndex(nn);
					nn++;
				}
			}
			newTracePatchPoints[nn] = tracePatchPoints[nTracePatchPoints - 1].resetIndex(nn);
			nn++;
			tracePatchPoints = (PatchPoint[]) StsMath.trimArray(newTracePatchPoints, nn);
			nTracePatchPoints = nn;
		}

		private final int intervalRound(float value, int min)
		{
			return min + Math.round(value - min);
		}

		private void initializePatchPointsList()
		{
			patchPointsList = new PatchPointLinkList(nTracePatchPoints);
		}

		/**
		 * Insert inactive row and col connections at the top and bottom of the patchPointsList.
		 * When scanning the patchPointsList, these have indices set to prevent going beyond the actual connections.
		 * @param prevColTrace trace in the same col, prev row as this trace
		 * @param prevRowTrace trace in the same row, prev col as this trace
		 */
		private void initializePatchPointsListConnections(TracePoints prevColTrace, TracePoints prevRowTrace)
		{
			if (prevRowTrace != null)
			{
				patchPointsList.first.rowConnection = new Connection(prevRowTrace.patchPointsList.first, patchPointsList.first);
				patchPointsList.last.rowConnection = new Connection(prevRowTrace.patchPointsList.last, patchPointsList.last);
			}
			if (prevColTrace != null)
			{
				patchPointsList.first.colConnection = new Connection(prevColTrace.patchPointsList.first, patchPointsList.first);
				patchPointsList.last.colConnection = new Connection(prevColTrace.patchPointsList.last, patchPointsList.last);
			}
		}

		private void reinitializeTraceIndices(TracePoints prevRowTrace, TracePoints prevColTrace)
		{
			if (trace != null) reinitializeTraceIndexing();
			if (prevRowTrace != null) prevRowTrace.reinitializeTraceIndexing();
			if (prevColTrace != null) prevColTrace.reinitializeTraceIndexing();
		}

		/**
		 * doubly linked list of PatchPoint[s].  PatchPoints in list are in increasing order of vertical index (slice).
		 * firstLink is at slice=-1, and lastLink is at slice=nPatchPoints (last real PatchPoint could only be at nPatchPoints-1).
		 */
		class PatchPointLinkList
		{
			/** first point in link list (connected to first actual point in list) */
			final PatchPoint first;
			/** last point in link list (connected to last actual point in list) */
			final PatchPoint last;
			/** last connected point in linked list just above current point */
			PatchPoint connectedPointAbove;
			/** connected point just below connectedPointAbove in linked list */
			PatchPoint connectedPointBelow;
			/** last connected point that was set; a convenient starting point for any search */
			PatchPoint currentConnectedPoint;

			PatchPointLinkList(int nTracePatchPoints)
			{
				first = new PatchPoint(row, col, -largeInt, -1.0f, 0.0f, POINT_ANY, -1);
				last = new PatchPoint(row, col, largeInt, seismicVolume.zMax + 1, 0.0f, POINT_ANY, nTracePatchPoints);
				first.next = last;
				last.prev = first;
				currentConnectedPoint = first;
				currentPoint = first;
				connectedPointAbove = first;
				connectedPointBelow = last;
				connectedPointAbove.next = last;
				connectedPointBelow.prev = first;
			}

			void reinitializeTraceIndexing()
			{
				connectedPointAbove = first;
				connectedPointBelow = first.next;
				currentConnectedPoint = first;
				currentPoint = first;
			}

			/**
			 * we have moved down to a new existing correlated patchPoint; set the interval to the one between this patchPoint and the point below
			 * @param connectedPoint the point at the top of the interval to be set
			 */
			void movePatchInterval(PatchPoint connectedPoint)
			{
				PatchPoint connectedPointClone = connectedPoint.clonedPoint;
				if(connectedPointClone != null)
					connectedPoint = connectedPointClone;
				connectedPointAbove = connectedPoint;
				connectedPointBelow = connectedPoint.next;
				currentPoint = connectedPoint;
				currentConnectedPoint = connectedPoint;
			}

			/**
			 * we are inserting this connectedPoint in an interval between connectedPointAbove and connectedPointBelow, either of which could be null
			 * meaning that it could be an openInterval with above and/or below undefined.  The interval (open or closed) is
			 * split into two subintervals and the current interval is set to the lower subinterval.
			 * @param connectedPoint point between pointAbove and pointBelow where interval is to be split into two subIntervals.
			 */
			void insert(PatchPoint connectedPoint)
			{
				if (connectedPoint == connectedPointAbove || connectedPoint == connectedPointBelow) return;
				connectedPointAbove.next = connectedPoint;
				connectedPoint.prev = connectedPointAbove;
				connectedPoint.next = connectedPointBelow;
				connectedPointBelow.prev = connectedPoint;

				connectedPointAbove = connectedPoint;
				connectedPointBelow = connectedPoint.next;
			}

			void debugPrintLinkList()
			{
				PatchPoint start = (PatchPoint) first;
				int n = 0;
				while (start != null)
				{
					System.out.println("patchPoint[" + n + "]" + start.toString());
					start = start.next;
					n++;
				}
				start = last;
				while (start != null)
				{
					System.out.println("patchPoint " + start.toString());
					start = start.prev;
					n--;
				}
			}
			/**
			 * Search up from startPoint for a point which has a row or col Connection.
			 * @param point point to start search from going up
			 * @param isRow point on this trace where we want to find/make a connection
			 * @return PatchPoint on the other row or col trace for the connection above the current search interval
			 */
			private PatchPoint getConnectedPatchPointAbove(PatchPoint point, boolean isRow)
			{
				int slice = point.slice;
				PatchPoint connectedPoint = currentConnectedPoint;
				try
				{
					if (connectedPoint.slice <= slice) // search down for connectedPoint until we bracket and return cpnnectedPoint just above point
					{
						while (connectedPoint.next.slice <= slice)
							connectedPoint = connectedPoint.next;
						return connectedPoint;
					}
					else // search up until we bracket and return connectedPoint just above point
					{
						while (connectedPoint.prev.slice > slice)
							connectedPoint = connectedPoint.prev;
						return connectedPoint;
					}
				}
				catch (Exception e)
				{
					StsException.outputWarningException(this, "getConnectedPatchPointAbove", e);
					return patchPointsList.first;
				}
			}

			/**
			 * rowConnections and colConnections have linked list of connections between this trace and the other row or col trace.
			 * For each list, connectionAbove and connectionBelow are the connections bounding our search interval on this trace.
			 * If isRow==true, return the patchPoint connected to the connectionBelow in rowConnections; otherwise
			 * return the patchPoint connected to the connectionBelow in colConnections.
			 * @param isRow point on this trace where we want to find/make a connection
			 * @return PatchPoint on the other row or col trace for the connection above the current search interval
			 */
			private PatchPoint getConnectedPatchPointBelow(PatchPoint point, boolean isRow)
			{
				int slice = point.slice;
				PatchPoint connectedPoint = currentConnectedPoint;
				try
				{
					if (connectedPoint.slice <= slice) // search down for connectedPoint until we bracket and return cpnnectedPoint just above point
					{
						while (connectedPoint.next.slice <= slice)
							connectedPoint = connectedPoint.next;
						return connectedPoint.next;
					}
					else // search up until we bracket and return connectedPoint just above point
					{
						while (connectedPoint.prev.slice > slice)
							connectedPoint = connectedPoint.prev;
						return connectedPoint;
					}
				}
				catch (Exception e)
				{
					StsException.outputWarningException(this, "getConnectedPatchPointAbove", e);
					return patchPointsList.last;
				}
			}
		}

		/**
		 * we have moved down to a new existing correlated patchPoint; set the interval to the one between this patchPoint and the point below
		 * @param connectedPoint the point at the top of the interval to be set
		 */
		void movePatchInterval(PatchPoint connectedPoint)
		{
			if (connectedPoint.hasConnection())
				patchPointsList.movePatchInterval(connectedPoint);
		}

		/**
		 * we are inserting this connectedPoint in an interval between connectedPointAbove and connectedPointBelow, either of which could be null
		 * meaning that it could be an openInterval with above and/or below undefined.  The interval (open or closed) is
		 * split into two subintervals and the current interval is set to the lower subinterval.
		 * @param connectedPoint point between pointAbove and pointBelow where interval is to be split into two subIntervals.
		 */
		void splitPatchInterval(PatchPoint connectedPoint)
		{
			if(connectedPoint.clonedPoint == null)
				patchPointsList.insert(connectedPoint);
		}

		void reinitializeTraceIndexing()
		{
			patchPointsList.reinitializeTraceIndexing();
		}

		/**
		 * currentPoint is the last currentPoint on this trace in the previous search operation, so is a good starting point for this search
		 * @param slice slice for which we want to find the nearest tracePoint
		 * @return the nearestTracePoint
		 */
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
					// if point is now below slice, then we have bracketed point: set and return currentPoint
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
				// didn't bracket, so slice is still below last point; return last point
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

					// if point is now above slice, then we have bracketed point: set and return currentPoint
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

		/**
		 * rowConnections and colConnections have linked list of connections between this trace and the other row or col trace.
		 * For each list, connectionAbove and connectionBelow are the connections bounding our search interval on this trace.
		 * If isRow==true, return the patchPoint connected to the connectionAbove in rowConnections; otherwise
		 * return the patchPoint connected to the connectionAbove in colConnections.
		 * @param isRow point on this trace where we want to find/make a connection
		 * @return PatchPoint on the other row or col trace for the connection above the current search interval
		 */
		private PatchPoint getOtherConnectedPatchPointAbove(boolean isRow)
		{

			PatchPoint connectionPointAbove = patchPointsList.connectedPointAbove;
			Connection connection = null;
			try
			{
				while (((connection = connectionPointAbove.getConnection(isRow)) == null) && connectionPointAbove.prev != null)
					connectionPointAbove = connectionPointAbove.prev;
				return connection.otherPoint;
			}
			catch (Exception e)
			{
				StsException.outputWarningException(this, "getOtherConnectedPatchPointAbove", e);
				return patchPointsList.first; // hack for now, don't allow exception!
			}
		}

		/**
		 * rowConnections and colConnections have linked list of connections between this trace and the other row or col trace.
		 * For each list, connectionAbove and connectionBelow are the connections bounding our search interval on this trace.
		 * If isRow==true, return the patchPoint connected to the connectionBelow in rowConnections; otherwise
		 * return the patchPoint connected to the connectionBelow in colConnections.
		 * @param isRow point on this trace where we want to find/make a connection
		 * @return PatchPoint on the other row or col trace for the connection above the current search interval
		 */
		private PatchPoint getOtherConnectedPatchPointBelow(boolean isRow)
		{
			PatchPoint connectionPointBelow = patchPointsList.connectedPointBelow;
			Connection connection = null;
			try
			{
				while (((connection = connectionPointBelow.getConnection(isRow)) == null) && connectionPointBelow.next != null)
					connectionPointBelow = connectionPointBelow.next;
				return connection.otherPoint;
			}
			catch (Exception e)
			{
				StsException.outputWarningException(this, "getOtherConnectedPatchPointBelow", e);
				return patchPointsList.last; // hack for now, don't allow exception!
			}
		}
	}

	private StsPatchGrid mergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint)
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

	private StsPatchGrid cantMergePatchGrids(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		StsPatchGrid otherPatchGrid = otherPatchPoint.getPatchGrid();
		StsPatchGrid newPatchGrid = newPatchPoint.getPatchGrid();
		// can't merge: find larger grid and add a clone of the overlapped point from the smaller grid to it

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


	private void checkAddPatchGridToRowGrids(StsPatchGrid patchGrid)
	{
		int patchID = patchGrid.id;
		if (patchGrid.rowGridAdded) return;
		boolean debug = patchGrid.debug();
		StsPatchGrid value = rowGrids.put(patchID, patchGrid); // if return is null, no value exists at this key
		patchGrid.rowGridAdded = true;
		if (debug)
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
		boolean debug = StsPatchGrid.debugPatchID != -1 && patchID == StsPatchGrid.debugPatchID;
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

	class PatchPoint implements Comparable<PatchPoint>, Cloneable
	{
		float value;
		float z = StsParameters.nullValue;
		byte pointType;
		StsPatchGrid patchGrid;
		int row;
		int col;
		int slice;

		PatchPoint next, prev;
		/** connection from this tracePoint to the tracePoint on the adjacent trace at this row, col-1 (same row) */
		Connection rowConnection;
		/** connection from this tracePoint to the tracePoint on the adjacent trace at this row-1, col (same col) */
		Connection colConnection;
		/** index of this point in the trace containing it */
		int traceIndex;
		/** correl factor between this point and next in row. Note that rowConnection is from this point back. */
		float rowCorrel;
		/** correl factor between this point and next in col.  Note that colConnection is from this point back. */
		float colCorrel;
		/** cloned point for debugging.  Point this point was cloned from if cloned. */
		PatchPoint clonedPoint;

		/** first point above which has a connected patch */
//		PatchPoint connectedPointAbove = null;

		/** first point below which has a connected patch */
//		PatchPoint connectedPointBelow = null;

		PatchPoint()
		{
		}

		/**
		 * constructor for first and last links in doubly-linked list of PatchPoints
		 */
		PatchPoint(int traceIndex)
		{
			this.traceIndex = traceIndex;
		}

		PatchPoint(int row, int col, int slice, float z, float value, byte pointType, int traceIndex)
		{
			this.row = row;
			this.col = col;
			this.slice = slice;
			this.z = z;
			this.value = value;
			this.pointType = pointType;
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

		/**
		 * A point needs to be cloned if it overlaps if by connected to it, the two newly connected grids overlap.
		 * So clone the point and clear connection data.  The single connection to this cloned point will be added.
		 * @return the cloned point
		 */
		protected PatchPoint clone()
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
				StsException.systemError(this, "clone");
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

		Integer hashCode(int nVolumeCols)
		{
			return new Integer(row * nVolumeCols + col);
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
			return col + row * nVolumeCols;
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
			return patchToString() + "r " + row + " c " + col + " s " + slice + " v " + value +
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

		public byte getPointType()
		{
			if (!useFalseTypes) return pointType;
			return StsTraceUtilities.coercedPointTypes[pointType];
		}

		public PatchPoint resetIndex(int index)
		{
			traceIndex = index;
			return this;
		}
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
		/** connected point on this trace */
		PatchPoint point;
		/** connected point on other trace */
		PatchPoint otherPoint;
		/** index in thisTrace.tracePoints array for this connection; -1 for topStart and nTracePoints for botEnd in link list */
		int thisTraceIndex;
		/** correlation between these two points; assigned to otherPoint location in either row or col direction */
		float correlation;

		Connection(int thisIndex)
		{
			thisTraceIndex = thisIndex;
		}

		Connection(PatchPoint newPoint, boolean isRow)
		{
			this.point = newPoint;
			otherPoint = newPoint.clone();
			if (isRow)
				otherPoint.col--;
			else
				otherPoint.row--;
		}

		Connection(PatchPoint otherPoint, PatchPoint newPoint, float correlation)
		{
			this(otherPoint, newPoint);
			this.correlation = correlation;
		}

		Connection(PatchPoint otherPoint, PatchPoint newPoint)
		{
			this.point = newPoint;
			this.otherPoint = otherPoint;
			thisTraceIndex = newPoint.traceIndex;
		}

		public String toString()
		{
			return " connection: other point " + otherPoint.toString() + " to " + point.toString() + " correl: " + correlation;
		}
	}

	/**
	 * main() used for debugPatchGrid test of TracePoints link list operations.
	 * @param args input nPoints followed by groups beginning with -1 and a series of integers indicating which points are randomly added.
	 * Each -1 value causes the insertion process to begin again at the top of the tracePoints.  Example:
	 * 14    -1  4  9     -1 1 12   -1 6
	 * 14 points with 3 passes where first 4 & 9 are added, then passes with 1 & 12, and finally 6 are added.
	 */
	public static void main(String[] args)
	{
		int nTotalPoints = Integer.parseInt(args[0]);
		int nValues = args.length - 1;
		int[] values = new int[nValues];
		int traceIndex = 0;
		for (int n = 0; n < nValues; n++)
			values[n] = Integer.parseInt(args[n + 1]);

		StsPatchVolume patchVolume = new StsPatchVolume();
		TracePoints tracePoints = patchVolume.testTraceLinkList(nTotalPoints);
		int nPass = 1;
		for (int n = 0; n < nValues; n++)
		{
			int index = values[n];
			if (index < 0) // start a new pass down the trace; reinitialize links
			{
				System.out.println("Pass: " + nPass++);
				tracePoints.reinitializeTraceIndexing();
				traceIndex = 0;
				continue;
			}
			// run down trace, checking at each patchPoint; if patchPoint has grid, adjust above and below indexes
			for (; traceIndex < index; traceIndex++)
			{
				// if this patchPoint is connected (has a grid), then reset trace above and below
				PatchPoint patchPoint = tracePoints.tracePatchPoints[traceIndex];
				if (patchPoint.getPatchGrid() != null)
				{
					tracePoints.movePatchInterval(patchPoint);
				}
			}
			// when we get to the index of the point being added, create fake grid for this point and set it
			PatchPoint patchPoint = tracePoints.tracePatchPoints[index];
			patchPoint.setPatchGrid(new StsPatchGrid());
			tracePoints.splitPatchInterval(patchPoint);
			traceIndex++;
		}
		// when completed, debugPatchGrid print the link list
		tracePoints.patchPointsList.debugPrintLinkList();
	}

	public TracePoints testTraceLinkList(int nPoints)
	{
		return new TracePoints(nPoints);
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
				// window size is even, so window ends with same point type as center; half-window size is windowSize/2.
				// we need to find windowSize/2 points above and below with same point type as window center (which is a max or min).
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
		if (StsPatchVolume.debug)
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
	 * then either delete it if it is a small point, or add it to volume set.
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
		if (debug)
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

				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawCol(gl, col, x, yMin, yInc, colorscale, false, displayCurvature);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID) ;
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
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID)
					gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
				patchGrid.drawRow(gl, row, y, xMin, xInc, colorscale, false, displayCurvature);
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID) ;
				gl.glLineWidth(getPatchVolumeClass().getEdgeWidth());
				// if (nFirst == -1) nFirst = n;
				if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID) break;
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
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawRow(gl, row, dirCoordinate, xMin, xInc, colorscale, true, displayCurvature);
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID) ;
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
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID)
						gl.glLineWidth(2 * getPatchVolumeClass().getEdgeWidth());
					patchGrid.drawCol(gl, col, dirCoordinate, yMin, yInc, colorscale, true, displayCurvature);
					if (drawPatchBold && patchGrid.id == StsPatchGrid.debugPatchInitialID) ;
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
					patchGrid.drawRow(gl, n, y, xMin, xInc, colorscale, is3d, displayCurvature);
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
		return patchPoint.row * nCols + patchPoint.col;
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