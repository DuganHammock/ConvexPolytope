(* ::Package:: *)

BeginPackage["DuganHammock`ConvexPolytope`"];

	(* Names *)
		ConvexPolytope;
		ConvexPolytopeChop;
		ConvexPolytopeData;
		ConvexPolytopeDimension;
		ConvexPolytopeDual;
		ConvexPolytopeEmbed;
		ConvexPolytopeFaceCounts;
		ConvexPolytopeFigure;
		ConvexPolytopeGraphics3D;
		ConvexPolytopeGraphics;
		ConvexPolytopeHalfSpaces;
		ConvexPolytopeIntersectingQ;
		ConvexPolytopeIntersection;
		ConvexPolytopeMeasure;
		ConvexPolytopeProject;
		ConvexPolytopeRandomPoint;
		ConvexPolytopeRegionMemberQ;
		ConvexPolytopeSlice;
		ConvexPolytopeTabView;
		ConvexPolytopeTransform;
		ConvexPolytopeTranslate;

		HalfSpacesAffineSpan;
		HalfSpacesBoundedQ;
		HalfSpacesCoordinateBounds;
		HalfSpacesFeasibleQ;
		HalfSpacesHyperPlanes;
		HalfSpacesInteriorPoint;
		HalfSpacesObjectiveBound;
		HalfSpacesObjectiveVertex;
		HalfSpacesRandomVertex;
		HalfSpacesReduce;
		HalfSpacesTransfer;

	(* usage *)
		ConvexPolytope::usage = "ConvexPolytope[ points ] returns the convex polytope of points, returned as a formatted data object ConvexPolytope[data].";
		ConvexPolytopeChop::usage = "ConvexPolytopeChop[polytope,normal,bound] returns the polytope as intersection with the half-space given by the normal and bound.";
		ConvexPolytopeData::usage = "ConvexPolytopeData[] creates an Association which defines the ConvexPolytope data format.";
		ConvexPolytopeDual::usage = "ConvexPolytopeDual[polytope] returns the dual of polytope via reciprocation about the unit sphere (centered at the polytope centroid).";
		ConvexPolytopeEmbed::usage = "ConvexPolytopeProject[polytope,subspace] returns the polytope as embedded into the subspace (assuming the polytope and subspace dimensions are equal). PolytopeEmbed[subspace] represents an operator that can be applied to a polytope.";
		ConvexPolytopeFaceCounts::usage = "PolytopeFaceCounts[polytope] outputs the tallies of face vertex counts for each face dimension.";
		ConvexPolytopeKeyMaster::usage = "ConvexPolytopeKeyMaster[polytope, method, (args\[Ellipsis])] is a handler for the ConvexPolytope data type which determines and calls the appropriate function as indicated by the given method.";
		ConvexPolytopeProject::usage = "ConvexPolytopeProject[polytope,subspace] returns the polytope as projected to the subspace. PolytopeProject[subspace] represents an operator that can be applied to a polytope.";
		ConvexPolytopeRandomPoint::usage = "ConvexPolytopeRandomPoint[ convexPolytope ] returns a random point in the interior of convexPolytope.\nConvexPolytopeRandomPoint[ convexPolytope, n ] returns a list of n random points in the interior of convexPolytope.";
		ConvexPolytopeSlice::usage = "ConvexPolytopeSlice[polytope,normal,bound] returns the hyper-planar slice of polytope.";
		ConvexPolytopeTransform::usage = "PolytopeTransform[polytope,matrix] returns the polytope with vertices and normals transformed by matrix.";
		ConvexPolytopeTranslate::usage = "PolytopeTranslate[polytope,vector] translates the polytope by vector, updating the polytope points and bounds.";

		ConvexPolytope::nokey = "\"`1`\" is not a valid key for data objects of type ConvexPolytope[\[Ellipsis]].";


	Begin["`Private`"];


(* ::Subsection:: *)
(*ConvexPolytope*)


	(* ConvexPolytope *)
		ConvexPolytopeData[] :=
			<|
				"Dimension" -> Minus[1],
				"EmbeddingDimension" -> Minus[1],
				"VertexCoordinates" -> {},
				"FaceVertexIndices" -> <||>,
				"SimplexVertexIndices" -> {},
				"AffineSpaceBasis" -> {},
				"NullSpaceNormals" -> {},
				"NullSpaceBounds" -> {},
				"HalfSpaceNormals" -> {},
				"HalfSpaceBounds" -> {}
			|>

		ConvexPolytope[] := ConvexPolytope[ConvexPolytopeData[]]
		ConvexPolytope[{}] := ConvexPolytope[]
		ConvexPolytope[ points_?MatrixQ ] :=
			Module[{
					affineDimension,
					convexPolytope
				},
				affineDimension = AffineSpanDimension[points];
				If[Equal[affineDimension,0],
					convexPolytope = ConvexPolytopeDimension0[points];
					Return[convexPolytope];
				];
				If[Equal[affineDimension,1],
					convexPolytope = ConvexPolytopeDimension1[points];
					Return[convexPolytope];
				];
				If[Greater[affineDimension,1],
					convexPolytope = ConvexPolytopeDimensionN[points];
					Return[convexPolytope];
				];
			]

		ConvexPolytopeDimension0[pointCoordinates_] := 
			Module[{
					embeddingDimension,
					vertex,
					data,
					convexPolytope
				},
				embeddingDimension = EmbeddingDimension[pointCoordinates];
				vertex = Part[pointCoordinates,1];
				data = ConvexPolytopeData[];
				AssociateTo[data,
					<|
						"Dimension" -> 0,
						"EmbeddingDimension" -> embeddingDimension,
						"VertexCoordinates" -> {vertex},
						"FaceVertexIndices" -> Association[Rule[0,{{1}}]],
						"SimplexVertexIndices" -> {{1}},
						"AffineSpaceBasis" -> {},
						"NullSpaceNormals" -> IdentityMatrix[embeddingDimension],
						"NullSpaceBounds" -> vertex,
						"HalfSpaceNormals" -> {},
						"HalfSpaceBounds" -> {}
					|>];
				Return[ConvexPolytope[data]];
			]

		ConvexPolytopeDimension1[pointCoordinates_] := 
			Module[{
					points,
					pointOrdering,
					embeddingDimension,
					vertex1,
					vertex2,
					normalVector,
					data,
					convexPolytope
				},
				points = pointCoordinates;
				embeddingDimension = EmbeddingDimension[points];
				(* lexicograph ordering gives extremal points *)
				pointOrdering = Ordering[N[points]];
				{vertex1, vertex2} = Part[points,Part[pointOrdering,{+1,-1}]];
				normalVector = NChop[Normalize[Subtract[vertex2,vertex1]]];
				data = ConvexPolytopeData[];
				AssociateTo[data,
					<|
						"Dimension" -> 1,
						"EmbeddingDimension" -> embeddingDimension,
						"VertexCoordinates" -> {vertex1,vertex2},
						"FaceVertexIndices" -> <|0->{{1},{2}},1->{{1,2}}|>,
						"SimplexVertexIndices" -> {{1,2}},
						"AffineSpaceBasis" -> {normalVector},
						"NullSpaceNormals" -> NChop[Orthogonalize[NullSpace[{normalVector}]]],
						"NullSpaceBounds" -> VectorProject[data["NullSpaceNormals"],vertex1] // NChop,
						"HalfSpaceNormals" -> {-normalVector,normalVector},
						"HalfSpaceBounds" -> {Dot[-normalVector,vertex1],Dot[normalVector,vertex2]}
					|>];
				convexPolytope = ConvexPolytope[data];
				Return[convexPolytope];
		]

		ConvexPolytopeDimensionN[pointsInput_] := 
			Module[{
					points,
					affineSpan,
					simplex,
					centroid,
					simplexes,
					vertexPointIndices,
					vertices,
					facets,
					normals,
					bounds,
					faces,
					normalSignatures,
					rulesPointToVertexIndices,
					data,
					convexPolytope
					},
					points = pointsInput;
					affineSpan = AffineSpan[points];
					{facets,normals,bounds,simplexes} = ConvexPolytopeSimplicialFacets[points];
					{vertices, facets, simplexes} = ConvexPolytopeReindexPointsToVertices[points, facets, simplexes];
					{faces, normals, bounds} = ConvexPolytopeMergeSimplicialFacets[facets,normals,bounds];
					{normals,bounds} = NChop[{normals,bounds}];
				(* output polytope *)
					data = ConvexPolytopeData[];
					data["Dimension"] = affineSpan["Dimension"];
					data["EmbeddingDimension"] = affineSpan["EmbeddingDimension"];
					data["VertexCoordinates"] = vertices;
					data["FaceVertexIndices"] = faces;
					data["SimplexVertexIndices"] = simplexes;
					data["AffineSpaceBasis"] = affineSpan["Basis"] // NChop;
					data["NullSpaceNormals"] = affineSpan["NullSpace"] // Select[Positive@*Chop@*Norm@*N] // NChop;
					data["NullSpaceBounds"] = VectorProject[vertices[[1]],data["NullSpaceNormals"]] // NChop;
					data["HalfSpaceNormals"] = normals // NChop;
					data["HalfSpaceBounds"] = bounds // NChop;
					convexPolytope = ConvexPolytope[data];
				(* for dim=2 polytopes *)
					If[data["Dimension"]==2, data["FaceVertexIndices",2] = {PolygonVertexOrdering[data["VertexCoordinates"]]};];
				Return[ConvexPolytope[data]];
			]
		ConvexPolytopeSimplicialFacets[pointCoordinatesInput_] :=
			Module[
				{
					pointCoordinates,
						pointIndices,
						embeddingDimension,
						centroid,
						affineSpan,
						affineSpaceDimension,
						affineNullSpaceDimension,
						affineSpaceBasis,
						affineNullSpaceBasis,
						affineSpaceBasePoint,
					simplexVertexPointIndices,
						simplexVertexCoordinates,
						lookupSimplexFacets,
						simplexes,
					vertexPointIndices,
						vertices,
					facets,
						normals,
						bounds,
					remPointIndices,
						remFacetPointHeights,
						remPointFacetHeights,
						remPointQs,
						remPointCoordinates,
					facetIndex,
						normal,
						bound,
						remHeights,
						maxHeight,
						maxPointIndices,
						newVertexPointIndex,
						newVertex,
						newVertexFacetSignatures,
						visibleFacetBs,
						visibleFacets,
						newSimplexes,
						newFacets,
						newNormals,
						newBounds
				},
				(* input *)
					pointCoordinates = pointCoordinatesInput;
					pointCoordinates = N[pointCoordinates];
					pointIndices = ListIndices[pointCoordinates];
				(* affine space *)
					affineSpan = AffineSpan[pointCoordinates];
					affineSpaceBasis = affineSpan["Basis"];

					affineSpaceDimension = affineSpan["Dimension"];
					affineNullSpaceBasis = affineSpan["NullSpace"];
					affineNullSpaceDimension = Length[affineNullSpaceBasis];
					affineSpaceBasePoint = affineSpan["Origin"];
				(* project to affine space *)
					pointCoordinates = VectorProject[pointCoordinates,affineSpaceBasis];

				(* initial affine simplex *)
					simplexVertexPointIndices = ConvexPolytopeExtremalSimplexVertexIndices[pointCoordinates];

					simplexVertexCoordinates = Part[pointCoordinates,simplexVertexPointIndices];
				(* recenter: initial simplex centroid \[RightTeeArrow] origin *)
					centroid = VectorCentroid[simplexVertexCoordinates];
					pointCoordinates = VectorTranslate[pointCoordinates,Minus[centroid]];
				(* simplex *)
					lookupSimplexFacets = Subsets[Range[Plus[1,affineSpaceDimension]],{affineSpaceDimension}];
				(* initial simplex *)
					simplexes = {simplexVertexPointIndices};
					vertexPointIndices = simplexVertexPointIndices;
					vertices = Part[pointCoordinates,vertexPointIndices];
				(* initial facets *)
					facets = Parts[simplexVertexPointIndices,lookupSimplexFacets];
					normals = VectorNormalize[Map[SimplexNormal,Parts[pointCoordinates,facets]]];
					bounds = NChop[MapThread[Dot,{normals,pointCoordinates[[facets[[All,1]]]]}]];
					(* orient normals *)
						normals *= Sign[bounds];
						bounds = Abs[bounds];
				(* build polytope simplices *)
					remPointIndices = Complement[pointIndices,simplexVertexPointIndices];
					While[Unequal[{},remPointIndices],
						(* remove interior remaining points *)
							If[remPointIndices=={},Break[];];
							remPointCoordinates = Part[pointCoordinates,remPointIndices];
							remFacetPointHeights = NChop[Subtract[Dot[normals,Transpose[remPointCoordinates]],bounds]];
							If[NoneTrue[Sign[remFacetPointHeights],EqualTo[1],2],Break[];];
							remPointFacetHeights = Transpose[remFacetPointHeights];
							remPointQs = Map[AnyTrue[Positive],Transpose[remFacetPointHeights]];
							remPointIndices = Pick[remPointIndices,remPointQs];
							remPointFacetHeights = Pick[remPointFacetHeights,remPointQs];
							remFacetPointHeights = Transpose[remPointFacetHeights];
						(* choose new vertex (extremal) *)
							facetIndex = First[Ordering[Map[Max,remFacetPointHeights],-1]];
							normal = normals[[facetIndex]];
							bound = bounds[[facetIndex]];
							remPointCoordinates = pointCoordinates[[remPointIndices]];
							remHeights = Subtract[Dot[remPointCoordinates,normal],bound];
							maxHeight = Max[remHeights];
							maxPointIndices = Pick[remPointIndices,remHeights,maxHeight];
							newVertexPointIndex = Part[maxPointIndices,First[Ordering[Part[pointCoordinates,maxPointIndices],Minus[1]]]];
							newVertex = Part[pointCoordinates,newVertexPointIndex];
						(* add new vertex to polytope *)
							vertexPointIndices = Append[vertexPointIndices,newVertexPointIndex];
							vertices = Append[vertices,newVertex];
						(* extrude visible facets to new simplexes *)
							visibleFacetBs = BooleNPositive[Subtract[Dot[normals,newVertex],bounds]];
							visibleFacets = Pick[facets,visibleFacetBs,1];
						(* keep non-visible {facets,normals,bounds} *)
							facets = Pick[facets,visibleFacetBs,0];
							normals = Pick[normals,visibleFacetBs,0];
							bounds = Pick[bounds,visibleFacetBs,0];
						(* get new facets from new simplexes - remove duplicated(interior) facets *)
							newSimplexes = ArrayPad[visibleFacets,{{0,0},{0,1}},newVertexPointIndex];
							newFacets = Flatten[Table[Parts[newSimplex,Rest[lookupSimplexFacets]],{newSimplex,newSimplexes}],1];
							newFacets = Sort[Map[Sort,newFacets]];
							newFacets = NonDuplicates[newFacets];
						(* add new simplexes to polytope *)
							newSimplexes = Map[Sort,newSimplexes];
							simplexes = Join[simplexes,newSimplexes];
						(* compute new affine normals, bounds *)
							newNormals = Map[SimplexNormal,Parts[pointCoordinates,newFacets]];
							newNormals = VectorNormalize[newNormals];
							newBounds = Dot[newNormals,newVertex];
							(* flip inward pointing normals *)
								newNormals *= Sign[newBounds];
								newBounds = Abs[newBounds];
						(* add new facets to polytope *)
							facets = Join[facets,newFacets];
							normals = Join[normals,newNormals];
							bounds = Join[bounds,newBounds];
						];
					(* recenter to (simplex centroid) \[Function] origin *)
						bounds = Plus[bounds,Dot[normals,centroid]];
					(* embed normals *)
						normals = VectorEmbed[normals,affineSpaceBasis];
					Return[{facets,normals,bounds,simplexes}];
				]

		ConvexPolytopeExtremalSimplexVertexIndices[ pointsInput_ ] := 
			Module[
				{
					points,
					embeddingDimension,
					vectors,
					pivot,
					nullSpaceBasis,
					extremalVertexIndex,
					extremalVectorIndex,
					extremalVector,
					simplexVertexIndices
				},
				(* input *)
					points = N[pointsInput];
					embeddingDimension = EmbeddingDimension[points];
				(* init *)
					nullSpaceBasis = IdentityMatrix[embeddingDimension];
					simplexVertexIndices = {};
				(* initial vertex *)
					pivot = First[points];
					extremalVertexIndex = First[Ordering[VectorNorm[VectorTranslate[points,Minus[pivot]]],Minus[1]]];
					pivot = Part[points,extremalVertexIndex];
					extremalVertexIndex = First[Ordering[VectorNorm[VectorTranslate[points,Minus[pivot]]],Minus[1]]];
					pivot = Part[points,extremalVertexIndex];
					vectors = VectorTranslate[points,Minus[pivot]];
				(* record initial vertex *)
					simplexVertexIndices = {extremalVertexIndex};
				Do[
					extremalVectorIndex = First[Ordering[VectorNormSquared[vectors],-1]];
					(* extremalVectorIndex = First[Ordering[vectors,-1]]; *)
					extremalVector = Part[vectors,extremalVectorIndex];
					If[Equal[Chop[Norm[extremalVector]],0],Break[];];
					simplexVertexIndices = Append[simplexVertexIndices,extremalVectorIndex];
					nullSpaceBasis = NullSpace[{extremalVector}];
					If[Equal[nullSpaceBasis,{}],Break[];];
					vectors = VectorProject[vectors,nullSpaceBasis];
				,{dimension,1,embeddingDimension}];
				Return[simplexVertexIndices];
			]

		ConvexPolytopeReindexPointsToVertices[pointCoordinates_, facetPointIndices_, simplexPointIndices_] := 
			Module[{
					vertexPointIndices,
					vertexIndices,
					vertexCoordinates,
					lookupPointToVertexIndices,
					facetVertexIndices,
					simplexVertexIndices
				},
				vertexPointIndices = Union[Flatten[facetPointIndices]];
				vertexIndices = ListIndices[vertexPointIndices];
				lookupPointToVertexIndices = ConstantArray[Null,Max[vertexPointIndices]];
				lookupPointToVertexIndices[[vertexPointIndices]] = vertexIndices;
				facetVertexIndices = Parts[lookupPointToVertexIndices,facetPointIndices];
				simplexVertexIndices = Parts[lookupPointToVertexIndices,simplexPointIndices];
				simplexVertexIndices = Sort[Map[Sort,simplexVertexIndices]];
				vertexCoordinates = Part[pointCoordinates,vertexPointIndices];
				Return[{vertexCoordinates, facetVertexIndices, simplexVertexIndices}];
			]

		ConvexPolytopeMergeSimplicialFacets[simplicialFacets_, simplicialFacetNormals_, simplicialFacetBounds_] := 
			Module[
				{
					polytopeDimension,
					facetComplexes,
					simplicialFacetIndexBins,
					mergedFacets,
					vertexFacetIncidences,
					GetIncidentMergedFacetIndices,
					faces
				},
				polytopeDimension = Length[First[simplicialFacets]];
				simplicialFacetIndexBins = GatherIndices[NRound[simplicialFacetNormals]];
				facetComplexes = Parts[simplicialFacets,simplicialFacetIndexBins];
				facetComplexes = Map[Sort,facetComplexes,2];
				facetComplexes = SortBy[facetComplexes,Flatten];
				mergedFacets = Apply[Union,facetComplexes,1];
				faces = Association[Table[Rule[n,{}],{n,0,polytopeDimension-1}]];
				faces[polytopeDimension-1] = mergedFacets;
				vertexFacetIncidences = IndexSubsetMemberships[mergedFacets];
				GetIncidentMergedFacetIndices = Function[face,Apply[Intersection,Part[vertexFacetIncidences,face]]];
				Module[
					{
						complexes,
						complexBoundarySimplexes,
						facetBoundarySimplexes,
						subSimplexes,
						subComplexes,
						subFaces
					},
					complexes = facetComplexes;
					Do[
						complexBoundarySimplexes = 
							Table[
								NonDuplicates[Map[Sort,Flatten[Subsets[Transpose[complexes[[complexIndex]]],{faceDimension}],{3,1}]]]
							,{complexIndex,Length[complexes]}];
						subSimplexes = Union[Flatten[complexBoundarySimplexes,1]];
						subComplexes = GatherBy[subSimplexes,GetIncidentMergedFacetIndices];
						subComplexes = SortBy[subComplexes,Flatten];
						subFaces = Apply[Union,subComplexes,1];
						faces[faceDimension-1] = subFaces;
						complexes = subComplexes;
					,{faceDimension,polytopeDimension-1,1,-1}];
				];
				Return[
					{
					faces,
					Part[simplicialFacetNormals,simplicialFacetIndexBins[[All,1]]],
					Part[simplicialFacetBounds,simplicialFacetIndexBins[[All,1]]]
					}
				];
			]
		ConvexPolytopeSimplicialFaceComplexes[ nFacesInput_ ] := 
			Module[{
					nFaces,
					nFaceComplexes,
					faces,
					face,
					boundaryFaceQs,
					boundaryComplexes,
					boundaryFaces,
					cornerVertexIndex,
					boundarySimplexes,
					faceSimplexes,
					faceComplexes,
					faceDimension
				},
				nFaces = nFacesInput;
				nFaceComplexes = nFaces;
				nFaceComplexes[0] = Map[List,nFaceComplexes[0],{-2}];
				nFaceComplexes[1] = Map[List,nFaceComplexes[1],{-2}];
				Do[
					nFaceComplexes[faceDimension] = 
					Table[
						face = nFaces[faceDimension][[faceIndex]];
						If[Equal[Length[face],Plus[1,faceDimension]],
							faceSimplexes = List[face];
							,
							cornerVertexIndex = Min[face];
							boundaryFaceQs = nFaces[faceDimension-1] // Map[Function[boundaryFace,And[Not[MemberQ[boundaryFace,cornerVertexIndex]],ContainsOnly[boundaryFace,face]]]];
							boundaryComplexes = Pick[nFaceComplexes[faceDimension-1],boundaryFaceQs];
							boundarySimplexes = Union[Apply[Join,boundaryComplexes]];
							faceSimplexes = ArrayPad[boundarySimplexes,{{0,0},{0,1}},cornerVertexIndex];
							faceSimplexes = Sort[Map[Sort,faceSimplexes]];
						];
						faceSimplexes
					,{faceIndex,Length[nFaces[faceDimension]]}];
				,{faceDimension,Complement[Keys[nFaces],{0,1}]}];
				Return[nFaceComplexes];
				]
		ConvexPolytopeHalfSpaces::NotFeasible = "halfspaces do not intersect";
		ConvexPolytopeHalfSpaces::NotBounded = "halfspace intersection not bounded";
		ConvexPolytopeHalfSpaces[{normalsInput_, boundsInput_}] := 
			Module[{
					embeddingDimension,
					halfSpaceNormals,
					halfSpaceBounds,
					hyperPlaneNormals,
					hyperPlaneBounds,
					halfSpaceNormalVectorSpace,
					halfSpaceNormalEmbeddedQ,
					feasibleQ,
					boundedQ,
					interiorPoint,
					affineSpan,
					recenterQ,
					affineEmbeddedQ,
					dualPoints,
					dualFacets,
					dualNormals,
					dualBounds,
					dualSimplexes,
					dualFacetIndices,
					vertices,
					convexPolytope
				},
				(* input *)
					embeddingDimension = EmbeddingDimension[normalsInput];
					halfSpaceNormals = NChop[N[normalsInput]];
					halfSpaceBounds = NChop[N[boundsInput]];
				
				(* recude half-spaces *)
					{halfSpaceNormals,halfSpaceBounds} = HalfSpacesReduce[{halfSpaceNormals,halfSpaceBounds}];

				(* project to span(halfspace normals) *)
					halfSpaceNormalVectorSpace = VectorSpan[halfSpaceNormals];
					halfSpaceNormalEmbeddedQ = Less[halfSpaceNormalVectorSpace["Dimension"],embeddingDimension];
					If[halfSpaceNormalEmbeddedQ,
						halfSpaceNormals = VectorProject[halfSpaceNormals,halfSpaceNormalVectorSpace["Basis"]];];
					
				(* check feasibility *)
					If[Not[HalfSpacesFeasibleQ[{halfSpaceNormals,halfSpaceBounds}]], 
						Message[ConvexPolytopeHalfSpaces::NotFeasible]; 
						Return[ConvexPolytope[]];
					];
				
				(* check bounded *)
					If[Not[HalfSpacesBoundedQ[{halfSpaceNormals,halfSpaceBounds}]],
						Message[ConvexPolytopeHalfSpaces::NotBounded]; 
						Return[ConvexPolytope[]];
					];
				
				(* affine span *)
					affineSpan = HalfSpacesAffineSpan[{halfSpaceNormals,halfSpaceBounds}];
				
				(* case: dimension \[Equal] 0 *)
					If[Equal[affineSpan["Dimension"],0], 
						Return[ConvexPolytope[{affineSpan["Origin"]}]];];
							
				(* recenter interiorPoint \[Function] Origin *)
					interiorPoint = HalfSpacesInteriorPoint[{halfSpaceNormals,halfSpaceBounds}];
					halfSpaceBounds = NChop[Subtract[halfSpaceBounds,Dot[halfSpaceNormals,interiorPoint]]];

				(* transfer halfspaces to affine space *)
					affineEmbeddedQ = Less[affineSpan["Dimension"],embeddingDimension];
					If[affineEmbeddedQ,
						{halfSpaceNormals,halfSpaceBounds} = HalfSpacesTransfer[{halfSpaceNormals,halfSpaceBounds},affineSpan["Basis"]];];
				(* get dual facets - spherical inversion *)
					dualPoints = Divide[halfSpaceNormals,halfSpaceBounds];

					{dualFacets,dualNormals,dualBounds,dualSimplexes} = ConvexPolytopeSimplicialFacets[dualPoints];
				(* select one simplex for each facet (after grouping) *)
					dualFacetIndices = Part[Values[PositionIndex[NRound[dualNormals]]],All,1];
					dualFacets = Part[dualFacets,dualFacetIndices];
				
				(* vertices *)
					vertices = Table[LinearSolve[Part[halfSpaceNormals,dualFacet],Part[halfSpaceBounds,dualFacet]],{dualFacet,dualFacets}];
				(* map vertices back into original vector/affine space *)
					If[affineEmbeddedQ, vertices = VectorEmbed[vertices,affineSpan["Basis"]];];
					vertices = VectorTranslate[vertices,interiorPoint];
					If[halfSpaceNormalEmbeddedQ, vertices = VectorEmbed[vertices,halfSpaceNormalVectorSpace["Basis"]]];
				(* chop numerical values *)
					vertices  = NChop[vertices];
					vertices = NumericalSort[vertices];
				(* polytope *)
					convexPolytope = ConvexPolytope[vertices];
				(* output *)
					Return[convexPolytope];
			]
		 
		ConvexPolytopeIntersection[a_,b_,c__] := ConvexPolytopeIntersection[ConvexPolytopeIntersection[a,b],c]
		ConvexPolytopeIntersection[ p1_ConvexPolytope, p2_ConvexPolytope ] := 
			Module[{aff1,aff2,aff,n,b,p},
				(* aff1 = AffineSpan[NChop[VectorCentroid[p1["VertexCoordinates"]]],p1["AffineSpaceBasis"]]; *)
				aff1 = AffineSpan[p1["VertexCoordinates"]];
				aff2 = AffineSpan[p2["VertexCoordinates"]];
				(* aff2 = AffineSpan[NChop[VectorCentroid[p2["VertexCoordinates"]]],p2["AffineSpaceBasis"]]; *)
				aff = AffineSpanIntersection[aff1,aff2];

				If[Equal[aff["Dimension"],Minus[1]],Return[ConvexPolytope[]];];
				If[Equal[aff["Dimension"],0],Return[ConvexPolytope[{aff["Origin"]}]];];
				n = Join[p1["HalfSpaceNormals"],p2["HalfSpaceNormals"]];
				b = Join[p1["HalfSpaceBounds"],p2["HalfSpaceBounds"]];
				(* p = ConvexPolytopeHalfSpaces[{n,b}]; *)
				(* re-center halfspaces at interior point *)
				(* b = Subtract[b,Dot[n,aff["Origin"]]]; *)
				{n,b} = HalfSpacesReduce[{n,b}];
				{n,b} = HalfSpacesTransfer[{n,b},aff["Basis"]];
				{n,b} = HalfSpacesReduce[{n,b}];
				n = NChop[n];
				b = NChop[b];
				p = ConvexPolytopeHalfSpaces[{n,b}];
				If[NonNegative[p["Dimension"]],
					p = ConvexPolytopeEmbed[p,aff["Basis"]];
					(* p = ConvexPolytopeTranslate[p,aff["Origin"]]; *)
				];
				Return[p];
			]

		ConvexPolytopeIntersectingQ[ p1_ConvexPolytope, p2_ConvexPolytope ] := 
			Module[{
					normals1, bounds1, signedBounds1,
					normals2, bounds2, signedBounds2,
					point1,
					point2,
					point1New,
					point2New,
					intersectingQ
				},
				normals1 = p1["HalfSpaceNormals"] // N;
				normals2 = p2["HalfSpaceNormals"] // N;
				bounds1 = p1["HalfSpaceBounds"] // N;
				bounds2 = p2["HalfSpaceBounds"] // N;
				signedBounds1 = Thread[{bounds1,-1}];
				signedBounds2 = Thread[{bounds2,-1}];
				
				point1 = p1["VertexCoordinates"][[1]] // N;
				point2 = p2["VertexCoordinates"][[1]] // N;

				intersectingQ = False;
				While[Not[intersectingQ],
					If[Or[VectorBoundedQ[point1, normals2, bounds2],VectorBoundedQ[point2, normals1, bounds1]],
						intersectingQ = True;
						Break[];
					];
					point1New = LinearProgramming[ Subtract[point2,point1], normals1, signedBounds1, Minus[Infinity]];
					point2New = LinearProgramming[ Subtract[point1,point2], normals2, signedBounds2, Minus[Infinity]];
					If[And[NEqual[point1,point1New],NEqual[point2,point2New]], 
						intersectingQ = False;
						Break[];
					];
				];
				Return[intersectingQ];
			];
		(* data *)
		ConvexPolytope[data_]["VertexCoordinates"|"Vertices"|"vertices","v"] := data["VertexCoordinates"]
		ConvexPolytope[data_]["FaceVertexIndices"|"FaceIndices"|"fvi"|"fi"] := data["FaceVertexIndices"]
		ConvexPolytope[data_]["FaceVertexCoordinates"|"FaceVertices"|"Faces"|"faces"|"fv"|"f",k_:All] := Parts[data["VertexCoordinates"],data["FaceVertexIndices",k]]
		ConvexPolytope[data_]["HalfSpaceNormals"|"FacetNormals"|"Normals"|"normals"|"n"|"hn"] := data["HalfSpaceNormals"]
		ConvexPolytope[data_]["HalfSpaceBounds"|"FacetBounds"|"Bounds"|"bounds"|"b"|"hb"] := data["HalfSpaceBounds"]
		ConvexPolytope[data_]["EmbeddingDimension"|"embdim"] := data["EmbeddingDimension"]
		ConvexPolytope[data_]["Dimension"|"dimension"|"dim"] := data["Dimension"]
		ConvexPolytope[data_]["Codimension"|"codimension"|"codim"] := Subtract[data["EmbeddingDimension"],data["Dimension"]]
		ConvexPolytope[data_]["AffineSpaceBasis"|"AffineSpace"|"AffineBasis"] := data["AffineBasis"]
		ConvexPolytope[data_]["NullSpaceNormals"|"NullSpaceBasis"|"NullSpace"|"NullVectors"] := data["NullSpaceNormals"]
		ConvexPolytope[data_]["NullSpaceBounds"|"NullBounds"] := data["NullSpaceBounds"]

		(* methods *)
		ConvexPolytope[data_]["Centroid"] := Mean[data["VertexCoordinates"]]
		ConvexPolytope[data_]["RandomPoint",n_] := ConvexPolytopeRandomPoint[ConvexPolytope[data],n]
		ConvexPolytope[data_]["Polygon"] := Polygon[ConvexPolytope[data]]
		ConvexPolytope[data_]["Dual"] := ConvexPolytopeDual[ConvexPolytope[data]]
		
		ConvexPolytope[data_][methodName_, args___] := ConvexPolytopeKeyMaster[ConvexPolytope[data],methodName,args]
		
		ConvexPolytope /: Part[ConvexPolytope[data_],partspec__] := Part[data,partspec]

		ConvexPolytopeKeyMaster[ConvexPolytope[data_], keyInput_, args___] := 
			Module[{key = keyInput},
				
				If[{"VertexCoordinates", "Vertices", "vertices", "v", "vc"} // MemberQ[key],
					Return[data["VertexCoordinates"]];];
				
				If[{"FaceVertexIndices", "fvi"} // MemberQ[key],
					Return[ConvexPolytopeFaces[ConvexPolytope[data],args]];];
				
				If[{"FaceVertexCoordiantes", "FaceVertices", "FaceCoordinates", "Faces", "faces", "f", "fv", "fvc", "fc"} // MemberQ[key],
					Return[Parts[data["VertexCoordinates"],ConvexPolytopeFaces[ConvexPolytope[data],args]]];];
				
				If[{"FacetVertexCoordinates", "FacetCoordinates", "FacetVertices", "Facets", "facets"} // MemberQ[key],
					Return[Parts[data["VertexCoordinates"],ConvexPolytopeFaces[ConvexPolytope[data],-1]]];];
					
				If[{"HalfSpaceNormals", "FacetNormals", "Normals", "normals", "n", "hn"} // MemberQ[key],
					Return[data["HalfSpaceNormals"]];];
				
				If[{"HalfSpaceBounds", "FacetBounds", "Bounds", "bounds", "b", "hb"} // MemberQ[key],
					Return[data["HalfSpaceBounds"]];];
				
				If[{"EmbeddingDimension", "AmbientDimension", "vdim", "embdim"} // MemberQ[key],
					Return[data["EmbeddingDimension"]];];
				
				If[{"Dimension", "AffineDimension", "dimension", "affdim", "dim"} // MemberQ[key],
					Return[data["Dimension"]];];

				If[{"Codimension", "codimension", "codim", "cd"} // MemberQ[key],
					Return[Subtract[data["EmbeddingDimension"],data["Dimension"]]];];
				
				If[{"FullDimensionQ", "DimensionFullQ"} // MemberQ[key],
					Return[Equal[data["Dimension"], data["EmbeddingDimension"]]];];
				
				If[{"AffineSpaceBasis", "AffineSpace", "affinespace", "AffineBasis", "affspace", "aff", "a"} // MemberQ[key],
					Return[data["AffineSpaceBasis"]];];
				
				If[{"AffineNullSpace", "affnull", "null", "perp", "NullSpace", "NullSpaceNormals", "NullSpaceBasis", "NullBasis"} // MemberQ[key],
					Return[data["NullSpaceNormals"]];];
				
				If[{"Centroid", "CentroidCoordinates", "centroid", "intp", "mean", "c"} // MemberQ[key],
					Return[Mean[data["VertexCoordinates"]]];];
				
				If[{"RandomPoint", "randp", "rand", "rp", "r"} // MemberQ[key],
					Return[ConvexPolytopeRandomPoint[ConvexPolytope[data], args]];];
				
				If[{"Translate", "translate"} // MemberQ[key],
					Return[ConvexPolytopeTranslate[ConvexPolytope[data], args]];];
				
				If[{"Embed", "embed"} // MemberQ[key],
					Return[ConvexPolytopeEmbed[ConvexPolytope[data], args]];];
					
				If[{"Project", "project", "Projection", "projection", "Proj", "proj"} // MemberQ[key],
					Return[ConvexPolytopeProject[ConvexPolytope[data], args]];];
					
				If[{"FaceCounts", "FaceTally", "f#", "fcount"} // MemberQ[key],
					Return[ConvexPolytopeFaceCounts[ConvexPolytope[data], args]];];
				
				If[{"Slice", "slice"} // MemberQ[key],
					Return[ConvexPolytopeSlice[ConvexPolytope[data], args]];];
				
				If[{"Figure", "fig"} // MemberQ[key],
					Return[ConvexPolytopeFigure[ConvexPolytope[data], args]];];
					
				If[{"TabView", "Tab", "tab"} // MemberQ[key],
					Return[ConvexPolytopeTabView[ConvexPolytope[data]]];];

				If[{"RegionMemberQ", "InteriorQ"} // MemberQ[key],
					Return[ConvexPolytopeRegionMemberQ[ConvexPolytope[data], args]];];
					
				If[{"Measure", "measure", "Volume", "volume", "vol"} // MemberQ[key],
					Return[ConvexPolytopeMeasure[ConvexPolytope[data], args]];];

				Message[ConvexPolytope::nokey, key];
				Return[Null];
			]

		ConvexPolytope /: MemberQ[point_][polytope_ConvexPolytope] := ConvexPolytopeRegionMemberQ[polytope,point]

		ConvexPolytopeFaces[ConvexPolytope[data_]] := data["FaceVertexIndices"]
		ConvexPolytopeFaces[ConvexPolytope[data_], k_Integer] := Lookup[data["FaceVertexIndices"],Mod[k,data["EmbeddingDimension"]]]

		ConvexPolytopeFaceCounts[ConvexPolytope[data_],faceDimensionsInput_:Null] := 
			Module[
				{
					faceDimensions,
					faceDimension,
					faceCodimension,
					polytopeDimension,
					faceVertexIndices,
					nfaces,
					tally,
					tallies,
					labels,
					label,
					result
				},
				faceDimensions = faceDimensionsInput;
				polytopeDimension = data["Dimension"];
				faceVertexIndices = data["FaceVertexIndices"];
				
				If[SameQ[faceDimensions,Null],
					faceDimensions = Range[0,Subtract[Length[faceVertexIndices],1]];
				];
				
				If[Not[ListQ[faceDimensions]],
					faceDimensions = List[faceDimensions];
				];
				
						
				tallies = 
				Table[
					tally = "\[EmptySet]";
					
					If[LessEqual[0,faceDimension,polytopeDimension],
						tally = MatrixForm[Tally[Map[Length,faceVertexIndices[[Plus[1,faceDimension]]]]]];
					];
					
					If[LessEqual[Minus[polytopeDimension],faceDimension,Minus[1]],
						faceCodimension = Abs[faceDimension];
						faceDimension = Subtract[polytopeDimension,faceCodimension];
						tally = MatrixForm[Tally[Map[Length,faceVertexIndices[[Plus[1,faceDimension]]]]]];
					];
					
					tally
					
				,{faceDimension,faceDimensions}];
				
				labels = 
					Table[
						label = "\[EmptySet]";
						
						If[LessEqual[0,faceDimension,polytopeDimension],
							label = ToString[faceDimension]<>"-faces";
						];
						
						If[LessEqual[Minus[polytopeDimension],faceDimension,Minus[1]],
							faceCodimension = Abs[faceDimension];
							faceDimension = Subtract[polytopeDimension,faceCodimension];
							label = ToString[faceDimension]<>"-faces";
						];
						
						label
				,{faceDimension,faceDimensions}];
				
				result = Grid[Transpose[MapThread[List,{labels,tallies}]], Dividers->1];
				
				Return[result];
			]

		ConvexPolytopeMeasure[ConvexPolytope[data_]] := 
			If[data["Dimension"]==0,1,
			If[Equal[{},data["SimplexVertexIndices"]],
				Undefined,
				Total[Map[SimplexMeasure,Parts[data["VertexCoordinates"],data["SimplexVertexIndices"]]]]
			]]

		ConvexPolytopeDimension[ConvexPolytope[data_]] := Lookup[data,"Dimension"]

		ConvexPolytopeFigure[ConvexPolytope[dataInput_],figureLabel_:"Polytope"] := 
			Module[{data,faceVertexIndices,faces,figure},
				data = dataInput;
				faceVertexIndices = data["FaceVertexIndices"];
				faces = Grid[Table[List["dim="<>ToString[faceDimension],Shallow[faceVertexIndices[[faceDimension+1]]]],{faceDimension,0,Length[faceVertexIndices]-1}],Rule[Dividers,Hue[0,0,0]],Rule[Alignment,{{Right,Left}}]];
				data = Map[Shallow,data];
				data["FaceVertexIndices"] = faces;
				figure = Grid[Apply[List,Normal[data],1],Rule[Dividers,1], Rule[Alignment,{{Right,Left}}]];
				Return[figure];
			]

		ConvexPolytopeOpenerView[ConvexPolytope[data_]] := 
			Module[{
					polytopeData,
					faceVertexIndices,
					facesTabView,
					dataOpenerView
				},
				polytopeData = data;
				dataOpenerView = polytopeData;
				dataOpenerView = MapAt[Shallow,dataOpenerView,{{"VertexCoordinates"},{"HalfSpaceNormals"},{"HalfSpaceBounds"},{"FaceVertexIndices",All}}];
				dataOpenerView["FaceVertexIndices"] = 
					OpenerView[{
						"FacesVertexIndices",
						Column[Table[
							OpenerView[{
								ToString[faceDimension]<>"-faces",
								dataOpenerView["FaceVertexIndices",faceDimension]
							}]
						,{faceDimension,0,polytopeData["Dimension"]-1}]]
					}];
				dataOpenerView = OpenerView[{"ConvexPolytope",Column[Map[OpenerView@*Apply[List],Normal[dataOpenerView]]]}];
				Return[dataOpenerView];
			]

		ConvexPolytopeTabView[ConvexPolytope[dataInput_]] := 
			Module[{data,faceVertexIndices,facesTabView,dataTabView},
				data = dataInput;
				faceVertexIndices = data["FaceVertexIndices"];
				facesTabView = 
					TabView[
						Table[Rule["dim="<>ToString[faceDimension]<>": ",Short[faceVertexIndices[[faceDimension+1]]]],{faceDimension,0,Length[faceVertexIndices]-1}]
					,ControlPlacement->Left];
				data = Map[Shallow,data];
				data["FaceVertexIndices"] = facesTabView;
				dataTabView = TabView[data,ControlPlacement->Left,ImageSize->512];
				Return[dataTabView];
			]

		ConvexPolytope /: TabView[ConvexPolytope[data_],options___] := TabView[data,options]
		ConvexPolytope /: MenuView[ConvexPolytope[data_],options___] := MenuView[data,options]
		ConvexPolytope /: Simplex[dilation_:1][polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Simplex, dilation]
		ConvexPolytope /: Simplex[polytope_ConvexPolytope, dilation_:1] := ConvexPolytopeGraphics3D[polytope, Simplex, dilation]
		ConvexPolytope /: Polygon[polytope_ConvexPolytope, dilation_:1] := ConvexPolytopeGraphics3D[polytope, Polygon, dilation]
		ConvexPolytope /: Polygon[dilation_:1][polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Polygon, dilation]
		ConvexPolytope /: Point[polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Point]
		ConvexPolytope /: Point[pointSize_][polytope_ConvexPolytope] := {AbsolutePointSize[pointSize],ConvexPolytopeGraphics3D[polytope,Point]}
		ConvexPolytope /: Line[polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Line]
		ConvexPolytope /: Line[thickness_][polytope_ConvexPolytope] := {AbsoluteThickness[thickness],ConvexPolytopeGraphics3D[polytope,Line]}
		ConvexPolytope /: Sphere[polytope_ConvexPolytope, radius_:0.1] := ConvexPolytopeGraphics3D[polytope, Sphere, radius]
		ConvexPolytope /: Sphere[radius_:0.1][polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Sphere, radius]
		ConvexPolytope /: Tube[polytope_ConvexPolytope, radius_:0.05, dilation_:1] := ConvexPolytopeGraphics3D[polytope, Tube, radius, dilation]
		ConvexPolytope /: Tube[radius_:0.05, dilation_:1][polytope_ConvexPolytope] := ConvexPolytopeGraphics3D[polytope, Tube, radius, dilation]

		ConvexPolytopeGraphics3D[polytope_ConvexPolytope?((#["Dimension"]==-1)&),___] := Null
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Point] := Point[VectorProject[polytope["VertexCoordinates"],3]]
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Sphere, radius_:1] := Sphere[VectorProject[polytope["VertexCoordinates"],3],radius]
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Line] := Line[Parts[VectorProject[polytope["VertexCoordinates"],3],polytope["FaceVertexIndices",1]]]
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Tube, radius_:.05, dilation_:1] := Tube[Map[VectorDilate[dilation],Parts[VectorProject[polytope["VertexCoordinates"],3],polytope["FaceVertexIndices",1]]],radius]
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Polygon, dilation_:1] := Map[Composition[Polygon,VectorDilate[dilation],PolygonVertexSort],Parts[VectorProject[polytope["VertexCoordinates"],3],polytope["FaceVertexIndices",2]]]
		ConvexPolytopeGraphics3D[polytope_ConvexPolytope, Simplex, dilation_:1] := 
			Module[{vertices,simplexes,graphics3D},
				If[Equal[-1,polytope["Dimension"]],Return[{}];];
				If[Equal[0,polytope["Dimension"]],Return[Simplex[{polytope["VertexCoordinates"][[1]]}]];];
				If[Equal[1,polytope["Dimension"]],Return[Simplex[{polytope["VertexCoordinates"]}]];];
				vertices = polytope["VertexCoordinates"] // VectorProject[3];
				simplexes = Parts[vertices3D,polytope["SimplexVertexIndices"]];
				If[Unequal[dilation,1],
					polygons3D = Map[VectorDilate[dilation],polygons3D];];
				graphics3D = Map[Simplex,polygons3D];
				Return[graphics3D];
			]

		ConvexPolytopeGraphics[ConvexPolytope[data_],Point] := 
			Module[{vertices2D,graphics2D},
				vertices2D = data["VertexCoordinates"] . IdentityMatrix[{data["EmbeddingDimension"],2}];
				graphics2D = Point[vertices2D];
				Return[graphics2D];
			]

		ConvexPolytopeGraphics[ConvexPolytope[data_],Disk,radius_:0.1] := 
			Module[{vertices2D,graphics2D},
				vertices2D = data["VertexCoordinates"] . IdentityMatrix[{data["EmbeddingDimension"],2}];
				graphics2D = Table[Disk[vertex,radius],{vertex,vertices2D}];
				Return[graphics2D];
			]

		ConvexPolytopeGraphics[ConvexPolytope[data_],Line] := 
			Module[{vertices3D,edges3D,graphics3D},
				vertices3D = data["VertexCoordinates"] . IdentityMatrix[{data["EmbeddingDimension"],2}];
				edges3D = Parts[vertices3D,data["FaceVertexIndices"][[2]]];
				graphics3D = Line[edges3D];
				Return[graphics3D];
			]

		ConvexPolytopeGraphics[ConvexPolytope[data_],Tube,width_:0.05] := 
			Module[{vertices,edges,graphics2D,normals,vectors,corners},
				vertices = data["VertexCoordinates"] . IdentityMatrix[{data["EmbeddingDimension"],2}];
				edges = Parts[vertices,data["FaceVertexIndices"][[2]]];
				vectors = Map[First@*Differences,edges];
				normals = VectorNormalize[Map[Cross,vectors]];
				corners = Table[edges[[index,1]] - normals[[index]]*width,{index,Length[edges]}];
				graphics2D = Table[Parallelogram[corners[[index]],{vectors[[index]],2*width*normals[[index]]}],{index,Length[edges]}];
				Return[graphics2D];
			]

		ConvexPolytopeGraphics[ConvexPolytope[data_],Polygon,dilation_:1] := 
			Module[{vertices2D,polygon2D,graphics2D},
				vertices2D = data["VertexCoordinates"] . IdentityMatrix[{data["EmbeddingDimension"],2}];
				polygon2D = PolygonOrient[vertices2D];
				If[Unequal[1,dilation],polygon2D = VectorDilate[polygon2D,dilation];];
				graphics2D = Polygon[polygon2D];
				Return[graphics2D];
			]

		ConvexPolytopeRandomPoint[polytope_ConvexPolytope, k_:1] :=
			Module[{
					vertices,

					randomWeights,
					randomPoints

				},
				polytope["VertexCoordinates"]
			]
		ConvexPolytopeRandomPoint[ConvexPolytope[data_],k_:1 ] := 
			Module[{ vertices, weights, randomPoints, randomPoint },
				vertices = N[data["VertexCoordinates"]];
				weights = RandomReal[{0,1},{arrayDimensions,Length[vertices]}]^data["Dimension"];
				weights /= weights . ConstantArray[1,Length[vertices]];
				randomPoints = Dot[weights,vertices];
				Return[randomPoints];
			]

		ConvexPolytopeRegionMemberQ[ConvexPolytope[data_],points_] := VectorBoundedQ[points,data["HalfSpaceNormals"],data["HalfSpaceBounds"]]
		ConvexPolytopeRegionMemberQ[ConvexPolytope[data_]][points_] := ConvexPolytopeRegionMemberQ[ConvexPolytope[data], points]

		ConvexPolytopeDual[polytope_] := 
			Module[{
					faceDimension,
					faceIndex,
					faceFacetIndices,
					dualPolytopeData,
					dualPolytope,
					polytopeRecenterQ,
					vertices,
					polytopeEmbedQ,
					dualBounds,
					pointFacetIncidences,
					origin,
					facetIndices,
					pointFacetAdj,
					pointFacetIndices,
					points,
					centroid,
					normals,
					bounds,
					facets,
					faces,
					affineDimension,
					embeddingDimension,
					dimension,
					affineSubspace,
					dualVertices,
					dualNormals,
					facetPoints,
					dualFaces
				},
				points = polytope["VertexCoordinates"];
				dimension = polytope["Dimension"];
				faces = polytope["FaceVertexIndices"];
				normals = polytope["HalfSpaceNormals"];
				bounds = polytope["HalfSpaceBounds"];
				points = points // N // Chop;
				normals = normals // N // Chop;
				bounds = bounds // N // Chop;
				facets = faces[dimension-1];
				centroid = Mean[points] // Chop;
				polytopeRecenterQ = AnyTrue[bounds,NonPositive];
				If[polytopeRecenterQ,
					points = VectorTranslate[points,Minus[centroid]];
					bounds = Subtract[bounds,Dot[normals,centroid]];
				];
				embeddingDimension = EmbeddingDimension[points];
				polytopeEmbedQ = Less[dimension,embeddingDimension];
				
				If[polytopeEmbedQ,
					affineSubspace = AffineSpan[points]["Basis"] // Orthogonalize // Chop;
					normals = VectorProject[normals,affineSubspace];
					points = VectorProject[points,affineSubspace];
				];
				
				facetIndices = Range[Length[facets]];
				vertices = points[[Flatten[faces[[1]]]]];
				dualVertices = Table[LinearSolve[points[[facet]],ConstantArray[1,Length[facet]]],{facet,facets}];
				pointFacetAdj = SparseArray[Rule[Flatten[MapIndexed[Composition[Tuples,List],facets],1],1]];
				pointFacetIndices = Map[Function[pointSubsetAdj,Pick[facetIndices,pointSubsetAdj,1]],pointFacetAdj];		
				pointFacetIndices = Round[pointFacetIndices];
				Monitor[
					dualFaces = Table[Apply[Intersection,pointFacetIndices[[faces[[faceDimension,faceIndex]]]]],{faceDimension,dimension},{faceIndex,Length[faces[[faceDimension]]]}];
					,Column[{ProgressIndicator[faceDimension,{1,dimension}],ProgressIndicator[faceIndex,{1,Length[faces[[faceDimension]]]}]}],1];
				dualFaces = Reverse[dualFaces];
				dualFaces = Association[Table[Rule[faceDimension,dualFaces[[faceDimension+1]]],{faceDimension,0,dimension-1}]];
				dualNormals = VectorNormalize[vertices];
				dualBounds = MapThread[Dot,{dualNormals,dualVertices[[dualFaces[[dimension,All,1]]]]}];
				
				If[polytopeEmbedQ,
					dualVertices = VectorEmbed[dualVertices,affineSubspace];
					dualNormals = VectorEmbed[dualNormals,affineSubspace];
					];
				
				If[polytopeRecenterQ,
					dualVertices = VectorTranslate[dualVertices,centroid];
					dualBounds += Dot[dualNormals,centroid];
					];

				dualVertices = NChop[dualVertices];
				dualNormals = NChop[dualNormals];
				dualBounds = NChop[dualBounds];
				dualPolytopeData = polytope // First;
				dualPolytopeData[["VertexCoordinates"]] = dualVertices;
				dualPolytopeData[["FaceVertexIndices"]] = dualFaces;
				dualPolytopeData[["HalfSpaceNormals"]] = dualNormals;
				dualPolytopeData[["HalfSpaceBounds"]] = dualBounds;
				dualPolytopeData[["Dimension"]] = dimension;
				dualPolytopeData[["SimplexVertexIndices"]] = {};
				dualPolytope = ConvexPolytope[dualPolytopeData];
				Return[dualPolytope];
			]

		ConvexPolytopeTranslate // ClearAll;
		ConvexPolytopeTranslate[ConvexPolytope[dataInput_], vector_] := 
			Module[{
					data,
					polytope,
					vertexCoordinates,
					pointCoordinates,
					halfSpaceNormals,
					halfSpaceBounds,
					NullSpaceNormals,
					NullSpaceBounds
				},
				data = dataInput;
				vertexCoordinates = data["VertexCoordinates"];
				pointCoordinates = data["PointCoordinates"];
				
				halfSpaceNormals = data["HalfSpaceNormals"];
				halfSpaceBounds = data["HalfSpaceBounds"];
				
				NullSpaceNormals = data["NullSpaceNormals"];
				NullSpaceBounds = data["NullSpaceBounds"];
				
				vertexCoordinates = Transpose[Plus[Transpose[vertexCoordinates],vector]];
				pointCoordinates = Transpose[Plus[Transpose[vertexCoordinates],vector]];
				halfSpaceBounds = Plus[halfSpaceBounds,Dot[halfSpaceNormals,vector]];
				
				If[Unequal[{},NullSpaceNormals],
					NullSpaceBounds = Plus[NullSpaceBounds,Dot[NullSpaceNormals,vector]];];
				
				data["VertexCoordinates"] = vertexCoordinates // NChop;
				data["HalfSpaceBounds"] = halfSpaceBounds // NChop;
				data["NullSpaceBounds"] = NullSpaceBounds // NChop;
				
				polytope = ConvexPolytope[data];
				Return[polytope];
			]
			
		ConvexPolytopeTranslate[vector_][polytope_] := ConvexPolytopeTranslate[polytope,vector]

		ConvexPolytopeEmbed[ConvexPolytope[dataInput_], basis_] := 
		Module[{
				data,
				embeddingDimension,
				basisNullSpace,
				vertexCoordinates,
				halfSpaceNormals,
				affineSpaceBasis,
				hyperPlaneNormals,
				hyperPlaneBounds,
				polytopeResult
			},
			data = dataInput;
			If[data["Dimension"]==-1, Return[ConvexPolytope[]];];
			
			embeddingDimension = EmbeddingDimension[basis];
			basisNullSpace = Orthogonalize[NullSpace[N[basis]]];
			
			vertexCoordinates = data["VertexCoordinates"];
			halfSpaceNormals = data["HalfSpaceNormals"];
			
			affineSpaceBasis = data["AffineSpaceBasis"];
			hyperPlaneNormals = data["NullSpaceNormals"];
			hyperPlaneBounds = data["NullSpaceBounds"];
			
			vertexCoordinates = VectorEmbed[vertexCoordinates,basis];
			If[Unequal[{},halfSpaceNormals], halfSpaceNormals = Dot[halfSpaceNormals,basis];];
			If[Unequal[{},affineSpaceBasis], affineSpaceBasis = Dot[affineSpaceBasis,basis];];
			If[Unequal[{},hyperPlaneNormals], hyperPlaneNormals = Dot[hyperPlaneNormals,basis];];
			
			If[Unequal[{},basisNullSpace],
				hyperPlaneNormals = Join[hyperPlaneNormals,basisNullSpace];
				hyperPlaneBounds = Join[hyperPlaneBounds,ConstantArray[0,Length[basisNullSpace]]];
			];
			
			data["VertexCoordinates"] = vertexCoordinates // NChop;
			data["HalfSpaceNormals"] = halfSpaceNormals // NChop;
			data["AffineSpaceBasis"] = affineSpaceBasis // NChop;
			data["NullSpaceNormals"] = hyperPlaneNormals // NChop;
			data["NullSpaceBounds"] = hyperPlaneBounds // NChop;
			data["EmbeddingDimension"] = embeddingDimension;
			
			polytopeResult = ConvexPolytope[data];
			
			Return[polytopeResult];
		]

		ConvexPolytopeEmbed[ basis_ ][ ConvexPolytope[data_] ] := ConvexPolytopeEmbed[ConvexPolytope[data],basis]

		ConvexPolytopeProject[ConvexPolytope[data_], basis_] := ConvexPolytope[NChop[VectorProject[data["VertexCoordinates"],basis]]]

		ConvexPolytopeTransform[basis_][polytope_ConvexPolytope] := ConvexPolytopeTransform[polytope, basis]
		ConvexPolytopeTransform[ConvexPolytope[polytopeDataInput_], basis_] := 
			Module[
				{
					polytope,
					embeddingDimension,
					basisNullSpace,
					vertexCoordinates,
					pointCoordinates,
					halfSpaceNormals,
					affineSpaceBasis,
					hyperPlaneNormals,
					polytopeResult
				},
				polytope = polytopeDataInput;
					
				vertexCoordinates = polytope["VertexCoordinates"];
				halfSpaceNormals = polytope["HalfSpaceNormals"];
				affineSpaceBasis = polytope["AffineSpaceBasis"];
				hyperPlaneNormals = polytope["NullSpaceNormals"];
				
				vertexCoordinates = VectorTransform[vertexCoordinates,basis];
				halfSpaceNormals = VectorTransform[halfSpaceNormals,basis];
				affineSpaceBasis = VectorTransform[affineSpaceBasis,basis];
				hyperPlaneNormals = VectorTransform[hyperPlaneNormals,basis];
				
				polytope["VertexCoordinates"] = vertexCoordinates;
				polytope["HalfSpaceNormals"] = halfSpaceNormals // NChop;
				polytope["AffineSpaceBasis"] = affineSpaceBasis // NChop;
				polytope["NullSpaceNormals"] = hyperPlaneNormals // NChop;
				
				polytopeResult = ConvexPolytope[polytope];
				
				Return[polytopeResult];
			]

		ConvexPolytopeEdgeSlice[{p1_,p2_},{v1_,v2_}] := 
			Switch[
				Sign[{v1,v2}],
				{-1,-1},{},
				{-1,0},{p2},
				{-1,1},{(p2*v1 - p1*v2)/(v1-v2)},
				{1,-1},{(p2*v1 - p1*v2)/(v1-v2)},
				{1,1},{},
				{1,0},{p2},
				{0,-1},{p1},
				{0,0},{p1,p2},
				{0,1},{p1}
			]

		ConvexPolytopeSlice[polytope_,normal_,value_] := 
			Module[{
					vertices,
					verticesN,
					vertexValues,
					vertexSigns,
					
					edgeVertexIndices,
					edgeIndices,
					edgeVertices,
					edgeVertexValues,
					edgeVertexSigns,
					
					slicedEdgeIndices,
					slicedEdgeVertexIndices,
					slicedEdgeVertices,
					slicedEdgeVertexValues,
					slicedEdgeVertexSigns,
					p1,p2,v1,v2,
					slicedVertices,
					sliceVertices,
					slicePolytope
				},
				vertices = polytope["vertices"];
				vertices = vertices // N;
				vertexValues = Chop[Subtract[Dot[vertices,normal],value]];
				vertexSigns = Sign[vertexValues];
				sliceVertices = Pick[vertices,vertexSigns,0];
				edgeVertexIndices = polytope["FaceVertexIndices"][1];
				edgeIndices = Range[Length[edgeVertexIndices]];
				edgeVertexSigns = Parts[vertexSigns,edgeVertexIndices];
				slicedEdgeIndices = Pick[edgeIndices,Map[Sort,edgeVertexSigns],{-1,1}];
				slicedEdgeVertexIndices = Part[edgeVertexIndices,slicedEdgeIndices];
					
				If[Positive[Length[slicedEdgeIndices]],
					slicedEdgeVertices = Parts[vertices,slicedEdgeVertexIndices];
					slicedEdgeVertexValues = Parts[vertexValues,slicedEdgeVertexIndices];
					{v1,v2} = Transpose[slicedEdgeVertexValues];
					{p1,p2} = Flatten[slicedEdgeVertices,{2}];
					slicedVertices = Divide[Chop[Subtract[Times[p2,v1],Times[p1,v2]]],Subtract[v1,v2]];
					,
					slicedVertices = {};
				];
				sliceVertices = Join[sliceVertices,slicedVertices];
				sliceVertices = sliceVertices // NDeleteDuplicates;
				slicePolytope = ConvexPolytope[sliceVertices];
				Return[slicePolytope];
			]
		ConvexPolytopeSlice[normal_,value_][polytope_] := ConvexPolytopeSlice[polytope,normal,value]

		ConvexPolytopeEdgeChop[{p1_,p2_},{v1_,v2_}] := 
			Switch[
				Sign[{v1,v2}],
				{-1,-1},{p1,p2},
				{-1,0},{p1,p2},
				{-1,1},{p1,(p2 v1-p1 v2)/(v1-v2)},
				{1,-1},{(p2 v1-p1 v2)/(v1-v2),p2},
				{1,1},Nothing,
				{1,0},{p2,p2},
				{0,-1},{p1,p2},
				{0,0},{p1,p2},
				{0,1},{p1,p1}
			]

		ConvexPolytopeChop[normals_,values_][polytope_] := ConvexPolytopeChop[polytope,normals,values]
		ConvexPolytopeChop[polytopeInput_,normalsInput_,boundsInput_] := 
			Module[{polytope,polytopeDimension,slicePolytope,polytopeNormals,polytopeNormalIndices,polytopeBounds,normals,normalIndices,bounds,edges,edgeValues,slicePoints,chopPolytope,normalIndex,normal,bound,points,newEdges,slicePointNormalIndices,slicePointIndices},
				polytope = polytopeInput;
				normals = normalsInput;
				normalIndices = ListIndices[normals];
				bounds = boundsInput;
				polytopeNormals = polytope["normals"];
				polytopeBounds = polytope["bounds"];
				polytopeDimension = polytope["dimension"];
				edges = Parts[polytope["vertices"],polytope["faces",1]];
				edges = edges // Map[NumericalSort];
				chopPolytope = ConvexPolytope[];
				Do[
					normal = normals[[normalIndex]];
					bound = bounds[[normalIndex]];
					edgeValues = Chop[Subtract[Dot[edges,normal],bound]];
					edges = MapThread[ConvexPolytopeEdgeChop,{edges,edgeValues}];
					If[normalIndex == Length[normals],Break[];];
					If[Equal[{},edges],Break[];];
					points = Flatten[edges,1] // NDeleteDuplicates;
					slicePoints = VectorMaximalByNormal[points,normal];			
					newEdges = Select[Table[VectorMaximalByNormal[slicePoints,N@normal],{normal,polytopeNormals}],EqualTo[2]@*Length];
					newEdges = newEdges // Map[NumericalSort];
					edges = Join[edges,newEdges] // NDeleteDuplicates;
				,{normalIndex,Length[normals]}];
				points = Flatten[edges,1] // NDeleteDuplicates;
				chopPolytope = ConvexPolytope[points];
				Return[chopPolytope];
			]

		HalfSpacesReduce[{normalsInput_?MatrixQ,boundsInput_?VectorQ}] := 
			Module[{
					normals,
					bounds
				},
				(* input *)
					normals = NChop[N[normalsInput]];
					bounds = NChop[N[boundsInput]];
					
				(* gather parallel normals *)
					Module[{indexBins},
						indexBins = GatherIndices[NRound[VectorNormalize[normals]]];
						normals = Part[normals,Part[indexBins,All,1]];
						bounds = Map[Min,Parts[bounds,indexBins]];
					];
				(* sort *)
					Module[{ordering},
						ordering = Ordering[normals];
						normals = Part[normals,ordering];
						bounds = Part[bounds,ordering];
					];

				Return[{normals,bounds}];
			]

		(* assume origin is interior *)
		HalfSpacesTransfer[{nInput_,bInput_},basisInput_] := 
			Module[{basis,n,b,nQs},
				(* input *)
					basis = N[basisInput];
					n = N[nInput];
					b = N[bInput];
				(* project normals to basis *)
					n = NChop[VectorProject[n,basis]];
				(* select (\[Eta],\[Beta])'s with non-trivial projection *)
					nQs = BooleNPositive[VectorNorm[n]];
					n = Pick[n,nQs,1];
					b = Pick[b,nQs,1];
				(* re-normalize *)
					b = Divide[b,VectorNorm[n]];
					n = VectorNormalize[n];
				Return[{n,b}];
			]

		HalfSpacesHyperPlanes[{halfSpaceNormalsInput_,halfSpaceBoundsInput_}] := 
			Module[{
					halfSpaceNormals,
					halfSpaceBounds,
					halfSpaceNormalBounds,
					hyperPlanesQ,
					hyperPlaneNormalBounds,
					indexBins,
					halfSpaceIndices,
					hyperPlaneIndices,
					hyperPlaneNormals,
					hyperPlaneBounds
				},
				halfSpaceNormals = halfSpaceNormalsInput;
				halfSpaceBounds = halfSpaceBoundsInput;
				hyperPlaneNormals = {};
				hyperPlaneBounds = {};
				
				(* check for antipodal half-spaces (hyperplanes) *)
					halfSpaceNormalBounds = Transpose[List[halfSpaceNormals,halfSpaceBounds]];
					indexBins = GatherIndicesBy[halfSpaceNormalBounds,Function[normalBound,Sort[NRound[List[-normalBound,normalBound]]]]];
					hyperPlanesQ = NoneTrue[indexBins,Composition[GreaterThan[1],Length]];
					If[Not[hyperPlanesQ],
						hyperPlaneNormals = {};
						hyperPlaneBounds = {};
						Return[{{halfSpaceNormals,halfSpaceBounds},{hyperPlaneNormals,hyperPlaneBounds}}];
					];
				
				halfSpaceIndices = Flatten[Select[indexBins,Composition[EqualTo[1],Length]]];
				hyperPlaneIndices = Part[Select[indexBins,Composition[EqualTo[2],Length]],All,1];
				hyperPlaneNormalBounds = halfSpaceNormalBounds[[hyperPlaneIndices]];
				halfSpaceNormalBounds = halfSpaceNormalBounds[[halfSpaceIndices]];

				hyperPlaneNormals = Part[hyperPlaneNormals,All,1];
				hyperPlaneBounds = Part[hyperPlaneNormalBounds,All,2];
				
				halfSpaceNormals = Part[halfSpaceNormals,All,1];
				halfSpaceBounds = Part[halfSpaceNormalBounds,All,2];
				
				Return[{{halfSpaceNormals,halfSpaceBounds},{hyperPlaneNormals,hyperPlaneBounds}}];
			]

		HalfSpacesObjectiveVertex[objectiveVector_, {halfSpaceNormals_, halfSpaceBounds_}] := 
			Module[{signedBounds,lowerBound,method,vertex},
				signedBounds = ArrayPad[Transpose[{halfSpaceBounds}],{{0,0},{0,1}},-1];
				method = "Simplex";
				lowerBound = DirectedInfinity[Minus[1]];
				vertex = Quiet[LinearProgramming[Minus[objectiveVector],halfSpaceNormals,signedBounds,lowerBound,Rule[Method,method]]];
				Return[vertex];
			]

		HalfSpacesObjectiveBound[objectiveVector_, {halfSpaceNormals_, halfSpaceBounds_}] := 
			Module[{vertex,bound},
				vertex = HalfSpacesObjectiveVertex[objectiveVector, {halfSpaceNormals, halfSpaceBounds}];
				bound = If[MemberQ[vertex,Indeterminate], Infinity, Dot[vertex,objectiveVector]];
				Return[bound];
			]

		HalfSpacesFeasibleQ[{normals_,bounds_}] := Not[SameQ[LinearProgramming,Head[Quiet[HalfSpacesRandomVertex[{normals,bounds}]]]]]

		HalfSpacesRandomVertex[{normals_,bounds_}] := HalfSpacesObjectiveVertex[RandomReal[{-1,1},EmbeddingDimension[normals]],{normals,bounds}]

		(* assume: normals are unit-length *)
		HalfSpacesInteriorPoint[{halfSpaceNormals_, halfSpaceBounds_}] := 
			Module[{
					embeddingDimension,
					signedCoordinateVectors,
					coordinateObjectiveVertices
				},
				embeddingDimension = EmbeddingDimension[halfSpaceNormals];
				signedCoordinateVectors = Join[IdentityMatrix[embeddingDimension],Minus[IdentityMatrix[embeddingDimension]]];
				vertices = Table[HalfSpacesObjectiveVertex[objectiveVector,{halfSpaceNormals,halfSpaceBounds}],{objectiveVector,signedCoordinateVectors}];
				interiorPoint = Mean[vertices];
				Return[interiorPoint];
			]

		HalfSpacesCoordinateBounds[{halfSpaceNormals_,halfSpaceBounds_}] := 
			Module[{embeddingDimension,coordinateIndex,coordinateIndices},
				embeddingDimension = EmbeddingDimension[halfSpaceNormals];
				Table[
						signValue*HalfSpacesObjectiveBound[Times[signValue,UnitVector[embeddingDimension,coordinateIndex]]
				,{halfSpaceNormals,halfSpaceBounds}]
				,{coordinateIndex,1,embeddingDimension}
				,{signValue,{-1,+1}}
				]
			]

		HalfSpacesAffineSpan[{normals_,bounds_}] := 
			Module[{
					halfSpaces,
					embeddingDimension,
					affineSpan,
					vertices,
					centroid,
					newVertices,
					dimension
				},
				embeddingDimension = EmbeddingDimension[normals];
				halfSpaces = NChop[N[{normals,bounds}]];
				vertices = {HalfSpacesRandomVertex[halfSpaces]};
				affineSpan = AffineSpan[vertices];
				Do[
					
					dimension = affineSpan["Dimension"];
					newVertices = NChop[Flatten[
						Table[
							HalfSpacesObjectiveVertex[Times[sign,objectiveVector],halfSpaces]
						,{sign,{-1,+1}}
						,{objectiveVector,affineSpan["NullSpace"]}
						],1]];
					vertices = NDeleteDuplicates[Join[vertices,newVertices]];
					centroid = Mean[N[vertices]];
					affineSpan = AffineSpan[Prepend[vertices,centroid]];
					If[Equal[dimension,affineSpan["Dimension"]],Break[];];
					If[Equal[embeddingDimension,affineSpan["Dimension"]],Break[];];
				,{embeddingDimension}];
				Return[affineSpan];
			]

		HalfSpacesBoundedQ[{halfSpaceNormals_,halfSpaceBounds_}] := 
			Module[{coordinateBounds,boundedQ},
				coordinateBounds = HalfSpacesCoordinateBounds[{halfSpaceNormals,halfSpaceBounds}];
				boundedQ = FreeQ[coordinateBounds,DirectedInfinity[1]];
				Return[boundedQ];
			];

		ConvexPolytope /: AffineSpan[ConvexPolytope[p_]] := AffineSpan[Mean[N[p["VertexCoordinates"]]],p["AffineSpaceBasis"]]
		ConvexPolytope /: VectorSpan[ConvexPolytope[p_]] := VectorSpan[p["AffineSpaceBasis"]]


(* ::Subsection::Closed:: *)
(*Boole*)


	(* Boole *)
		(* utility functions for Boolean operations on binary arrays *)
		(* Boole => Boolean *)
			Boolean[b_] := Positive[b]
		(* Boole | Logic *)
			(* Boole | Logic | Not *)
				BooleNot[b_] := Subtract[1,b]
			(* Boole | Logic | Compare *)
				BooleAnd[x_,y__] := BitAnd[x,y]
				BooleNand[x_,y__] := BooleNot[BooleAnd[x,y]]
				BooleNor[x_,y__] := BooleNot[BooleOr[x,y]]
				BooleOr[x_,y__] := BitOr[x,y]
				BooleXor[x_,y__] := BitXor[x,y]
			(* Boole | Logic | Test *)
				AllBoole[b__] := Equal[1,Min[b]]
				AnyBoole[b__] := Equal[1,Max[b]]
				NoneBoole[b__] := Equal[0,Max[b]]
		(* Boole | Sign *)
			BoolePositive[x_] := Ramp[Sign[x]]
			BooleNonPositive[x_] := UnitStep[Minus[x]]
			BooleEqualToZero[x_] := BooleNot[BooleUnequalToZero[x]]
			BooleUnequalToZero[x_] := Unitize[x]
			BooleNegative[x_] := BooleNot[BooleNonNegative[x]]
			BooleNonNegative[x_] := UnitStep[x]
		(* Boole | Compare *)
			BooleEqual[x_,y_] := BooleEqualToZero[Subtract[x,y]]
			BooleGreater[x_,y_] := BoolePositive[Subtract[x,y]]
			BooleGreaterEqual[x_,y_] := BooleNonNegative[Subtract[x,y]]
			BooleLess[x_,y_] := BooleNegative[Subtract[x,y]]
			BooleLessEqual[x_,y_] := BooleNonPositive[Subtract[x,y]]
			BooleUnequal[x_,y_] := BooleUnequalToZero[Subtract[x,y]]
			(* operators *)
				BooleEqualTo[y_][x_] := BooleEqual[x,y]
				BooleGreaterThan[y_][x_] := BooleGreater[x,y]
				BooleGreaterEqualThan[y_][x_] := BooleGreaterEqual[x,y]
				BooleLessThan[y_][x_] := BooleLess[x,y]
				BooleLessEqualThan[y_][x_] := BooleLessEqual[x,y]
				BooleUnequalTo[y_][x_] := BooleUnequal[x,y]


(* ::Subsection::Closed:: *)
(*List*)


	(* List *)
		(* utility functions for lists *)
		(* Parts *)
			Parts[list_, i_Integer] := Part[list,i]
			Parts[list_, {indices___Integer}] := Part[list,{indices}]
			Parts[list_, indexLists_List] := Map[Function[Parts[list,#]], indexLists]
			Parts[list_][indices_] := Parts[list,indices]
		(* Indices *)
			ListIndices[list_List] := Range[Length[list]]
			ElementIndices[x_,element_] := Pick[ListIndices[x],x,element]
			SublistIndices[list_,sublist_] := Subtract[Part[Take[Values[PositionIndex[Join[sublist,list]]],Length[sublist]],All,2],Length[sublist]]
		(* IndexSubset *)
			IndexSubsetAdjacencyMatrix[indexSubsets_, maxIndex_:0] := Normal[IndexSubsetAdjacencyMatrixSparse[indexSubsets, maxIndex]]
			IndexSubsetAdjacencyMatrixSparse[indexSubsets_, maxIndex_:0] := SparseArray[Rule[Flatten[MapIndexed[Composition[Tuples,List],indexSubsets],1],1],{Max[maxIndex,indexSubsets],Length[indexSubsets]}]
			IndexSubsetMemberships[pointIndexSubsets_] := Map[Function[pointSubsetAdj,Pick[ListIndices[pointIndexSubsets],pointSubsetAdj,1]],IndexSubsetAdjacencyMatrixSparse[pointIndexSubsets]]
		(* Gather *)
			GatherIndices[x_] := Values[PositionIndex[x]]
			GatherIndicesBy[x_, f_] := Values[PositionIndex[Map[f,x]]]
			GatherByList[x_, y_] := Parts[x,GatherIndices[y]]
		(* Duplicates *)
			(* Duplicates | Indices *)
				DuplicateIndexSets[x_] := Cases[PositionIndex[x],{i_,j__} :> {i,j}]
				DuplicateIndices[x_] := Apply[Union,DuplicateIndexSets[x]]
				FirstDuplicateIndices[x_] := Cases[PositionIndex[x],{i_,j__} :> i]
				NonDuplicateIndices[x_] := Cases[PositionIndex[x],{i_}:>i]
				DeleteDuplicateIndices[x_] := Part[Values[PositionIndex[x]],All,1]
			(* Duplicates *)
				Duplicates[x_] := Cases[Tally[x],RuleDelayed[{e_,Except[1]},e]]
				NonDuplicates[x_] := Cases[Tally[x],RuleDelayed[{e_,1},e]]
				FirstDuplicates[x_] := Part[x,FirstDuplicateIndices[x]]
				DuplicateSets[x_] := Cases[Tally[x],{i_}:>i]
				(* DeleteDuplicates *)
		(* Set *)
			(* Set | Empty *)
				EmptySet[] := {}
				EmptySetQ[x_] := Equal[0,Length[x]]
				SetNonEmptyQ[x_] := Positive[Length[x]]
			(* Set | Indices *)
				IntersectionIndices[x_,y_] := SublistIndices[x,Intersection[x,y]]
				ComplementIndices[x_,y_] := SublistIndices[x,Complement[x,y]]
		(* Maximal *)
			(* Maximal | Indices *)
				MaximalIndex[x_] := First[Ordering[x,-1]]
				MaximalIndices[x_] := Pick[ListIndices[x], x, MaximalElement[x]]
			(* Maximal | Elements *)
				MaximalElement[x_] := Part[x, MaximalIndex[x]]
				MaximalElements[x_] := Part[x, MaximalIndices[x]]
		(* Minimal *)
			(* Minimal | Indices *)
				MinimalIndex[x_] := First[Ordering[x,1]]
				MinimalIndices[x_] := Pick[ListIndices[x],x,MinimalElement[x]]
			(* Minimal | Elements *)
				MinimalElement[x_] := Part[x,MinimalIndex[x]]
				MinimalElements[x_] := Part[x,MinimalIndices[x]]
		(* Sign *)
			(* All | Sign *)
				AllNegative[x_] := AllBoole[BooleNegative[x]]
				AllNonNegative[x_] := AllBoole[BooleNonNegative[x]]
				AllUnequalToZero[x_] := AllBoole[BooleUnequalToZero[x]]
				AllPositive[x_] := AllBoole[BoolePositive[x]]
				AllEqualToZero[x_] := AllBoole[BooleEqualToZero[x]]
			(* Any | Sign *)
				AnyNegative[x_] := AnyBoole[BooleNegative[x]]
				AnyNonNegative[x_] := AnyBoole[BooleNonNegative[x]]
				AnyUnequalToZero[x_] := AnyBoole[BooleUnequalToZero[x]]
				AnyPositive[x_] := AnyBoole[BoolePositive[x]]
				AnyEqualToZero[x_] := AnyBoole[BooleEqualToZero[x]]
			(* None | Sign *)
				NoneNegative[x_] := NoneBoole[BooleNegative[x]]
				NoneNonNegative[x_] := NoneBoole[BooleNonNegative[x]]
				NoneUnequalToZero[x_] := NoneBoole[BooleUnequalToZero[x]]
				NonePositive[x_] := NoneBoole[BoolePositive[x]]
				NoneEqualToZero[x_] := NoneBoole[BooleEqualToZero[x]]
		(* Compare *)
			(* for comparing numerical lists *)
			(* List | Compare *)
				ListEqual[x_,y_] := Boolean[BooleEqualToZero[Subtract[x,y]]]
				ListGreater[x_,y_] := Boolean[BooleGreater[x,y]]
				ListGreaterEqual[x_,y_] := Boolean[BooleGreaterEqual[x,y]]
				ListLess[x_,y_] := Boolean[BooleLess[x,y]]
				ListLessEqual[x_,y_] := Boolean[BooleLessEqual[x,y]]
				ListUnequal[x_,y_] := Boolean[BooleUnequal[x,y]]
				(* operators *)
					ListEqualTo[y_][x_] := ListEqual[x,y]
					ListGreaterThan[y_][x_] := ListGreater[x,y]
					ListGreaterEqualThan[y_][x_] := ListGreaterEqual[x,y]
					ListLessThan[y_][x_] := ListLess[x,y]
					ListLessEqualThan[y_][x_] := ListLessEqual[x,y]
					ListUnequalTo[y_][x_] := ListUnequal[x,y]
			(* All | Compare *)
				AllEqual[x_,y_] := AllBoole[BooleEqual[x,y]]
				AllGreater[x_,y_] := AllBoole[BooleGreater[x,y]]
				AllGreaterEqual[x_,y_] := AllBoole[BooleGreaterEqual[x,y]]
				AllLess[x_,y_] := AllBoole[BooleLess[x,y]]
				AllLessEqual[x_,y_] := AllBoole[BooleLessEqual[x,y]]
				AllUnequal[x_,y_] := AllBoole[BooleUnequal[x,y]]
				(* operators *)
					AllEqualTo[y_][x_] := AllEqual[x,y]
					AllGreaterThan[y_][x_] := AllGreater[x,y]
					AllGreaterEqualThan[y_][x_] := AllGreaterEqual[x,y]
					AllLessThan[y_][x_] := AllLess[x,y]
					AllLessEqualThan[y_][x_] := AllLessEqual[x,y]
					AllUnequalTo[y_][x_] := AllUnequal[x,y]
			(* Any | Compare *)
				AnyEqual[x_,y_] := AnyBoole[BooleEqual[x,y]]
				AnyGreater[x_,y_] := AnyBoole[BooleGreater[x,y]]
				AnyGreaterEqual[x_,y_] := AnyBoole[BooleGreaterEqual[x,y]]
				AnyLess[x_,y_] := AnyBoole[BooleLess[x,y]]
				AnyLessEqual[x_,y_] := AnyBoole[BooleLessEqual[x,y]]
				AnyUnequal[x_,y_] := AnyBoole[BooleUnequal[x,y]]
				(* operators *)
					AnyEqualTo[y_][x_] := AnyEqual[x,y]
					AnyGreaterThan[y_][x_] := AnyGreater[x,y]
					AnyGreaterEqualThan[y_][x_] := AnyGreaterEqual[x,y]
					AnyLessThan[y_][x_] := AnyLess[x,y]
					AnyLessEqualThan[y_][x_] := AnyLessEqual[x,y]
					AnyUnequalTo[y_][x_] := AnyUnequal[x,y]
			(* None | Compare *)
				NoneEqual[x_,y_] := NoneBoole[BooleEqual[x,y]]
				NoneGreater[x_,y_] := NoneBoole[BooleGreater[x,y]]
				NoneGreaterEqual[x_,y_] := NoneBoole[BooleGreaterEqual[x,y]]
				NoneLess[x_,y_] := NoneBoole[BooleLess[x,y]]
				NoneLessEqual[x_,y_] := NoneBoole[BooleLessEqual[x,y]]
				NoneUnequal[x_,y_] := NoneBoole[BooleUnequal[x,y]]
				(* operators *)
					NoneEqualTo[y_][x_] := NoneEqual[x,y]
					NoneGreaterThan[y_][x_] := NoneGreater[x,y]
					NoneGreaterEqualThan[y_][x_] := NoneGreaterEqual[x,y]
					NoneLessThan[y_][x_] := NoneLess[x,y]
					NoneLessEqualThan[y_][x_] := NoneLessEqual[x,y]
					NoneUnequalTo[y_][x_] := NoneUnequal[x,y]


(* ::Subsection::Closed:: *)
(*Numeric*)


	(* Numeric *)
		(* N | Precision *)
			$NTolerance = 10.^-10;
		(* N | Chop *)
			NChop[{}] := {}
			NChop[x_?ArrayQ] := Threshold[x,{"Hard",$NTolerance}]
			NChop[x_List] := Map[NChop,x]
			NChop[x_] := N[Chop[x,$NTolerance]]
		(* N | Pick *)
			NPick[list_, sel_, pattern_] := Pick[list, NRound[sel], NRound[pattern]]
		(* N | Round *)
			NRound[x_] := Round[x,$NTolerance]
		(* Unit *)
			NonUnitize[x_] := Subtract[1,Unitize[x]]
			UnitDelta[x_] := NonUnitize[x]
			InvertUnits[x_] := With[{xUnital = Unitize[x]},Times[Power[Plus[x,Subtract[1,xUnital]],-1],xUnital]]
		(* N | Unit *)
			NUnitize[x_] := Unitize[NChop[x]]
			NNonUnitize[x_] := Subtract[1,NUnitize[x]]
			NUnitDelta[x_] := NNonUnitize[x]
			NUnitStep[x_] := UnitStep[Plus[x,$NTolerance]]
			NInvertUnits[x_] := With[{xNUnital = NUnitize[x]},Times[Power[Plus[x,Subtract[1,xNUnital]],-1],xNUnital]]
		(* N | Boole *)
			(* N | Boole | Sign *)
				BooleNEqualToZero[x_] := NUnitDelta[x]
				BooleNNegative[x_] := UnitStep[Subtract[Minus[$NTolerance],x]]
				BooleNNonNegative[x_] := UnitStep[Plus[x,$NTolerance]]
				BooleNNonPositive[x_] := UnitStep[Subtract[$NTolerance,x]]
				BooleNPositive[x_] := UnitStep[Plus[x,Minus[$NTolerance]]]
				BooleNUnequalToZero[x_] := NUnitize[x]
			(* N | Boole | Compare *)
				BooleNEqual[x_,y_] := BooleNEqualToZero[Subtract[x,y]]
				BooleNGreater[x_,y_] := BooleNPositive[Subtract[x,y]]
				BooleNGreaterEqual[x_,y_] := BooleNNonNegative[Subtract[x,y]]
				BooleNLess[x_,y_] := BooleNNegative[Subtract[x,y]]
				BooleNLessEqual[x_,y_] := BooleNNonPositive[Subtract[x,y]]
				BooleNUnequal[x_,y_] := BooleNEqualToZero[Subtract[x,y]]
		(* N | Gather *)
			NGather[x_] := GatherByList[x,NRound[x]]
			NGatherByList[x_,y_] := GatherByList[x,NRound[y]]
			NGatherIndices[x_] := GatherIndices[NRound[x]]
		(* N | Duplicates *)
			(* N | Duplicates | Indices *)
				NDeleteDuplicateIndices[x_] := DeleteDuplicateIndices[NRound[x]]
				NDuplicateIndexSets[x_] := DuplicateIndexSets[NRound[x]]
				NDuplicateIndices[x_] := DuplicateIndices[NRound[x]]
				NFirstDuplicateIndices[x_] := FirstDuplicateIndices[NRound[x]]
				NNonDuplicateIndices[x_] := NonDuplicateIndices[NRound[x]]
			(* N | Duplicates *)
				NDeleteDuplicates[x_] := Part[x,NDeleteDuplicateIndices[x]]
				NDuplicates[x_] := Part[x,NDuplicateIndices[x]]
				NDuplicateSets[x_] := Part[x,NDuplicateIndexSets[x]]
				NFirstDuplicates[x_] := Part[x,NFirstDuplicateIndices[x]]
				NNonDuplicates[x_] := Part[x,NNonDuplicateIndices[x]]
			(* N | Union *)
				NUnion[x_] := NumericalSort[NDeleteDuplicates[x]]
		(* N | Set *)
			(* N | Complement *)
				NComplementIndices[x_,y_] := ComplementIndices[NRound[x],NRound[y]]
				NComplement[x_,y_] := Part[x,NComplementIndices[x,y]]
				NComplement[x_,y_,z__] := NComplement[NComplement[x,y],z]
				NMemberQ[list_, form_] := MemberQ[NRound[list],NRound[form]]
			(* N | Intersection *)
				NIntersectionIndices[x_,y_] := IntersectionIndices[NRound[x],NRound[y]]
				NIntersection[x_,y_] := Part[x,NIntersectionIndices[x,y]]
				NIntersection[x_,y_,z__] := NIntersection[NIntersection[x,y],z]
		(* N | Maximal *)
			NMaximalElement[x_] := NRound[Max[x]]
			NMaximalElements[x_] := Part[x,MaximalIndices[x]]
			NMaximalIndex[x_] := First[Ordering[x,-1]]
			NMaximalIndices[x_] := Pick[ListIndices[x], x, MaximalElement[x]]
		(* N | Minimal *)
			NMinimalElement[x_] := Part[x,MinimalIndex[x]]
			NMinimalElements[x_] := Part[x,MinimalIndices[x]]
			NMinimalIndex[x_] := First[Ordering[x,-1]]
			NMinimalIndices[x_] := Pick[ListIndices[x],x,MinimalElement[x]]
		(* N | Sign *)
			(* N | Sign *)
				NSign[x_] := Sign[NChop[x]]
			(* N | Sign | Test *)
				NPositive[x_] := Boolean[BooleNPositive[x]]
				NNonPositive[x_] := Boolean[BooleNNonPositive[x]]
				NEqualToZero[x_] := Boolean[BooleNUnequalToZero[x]]
				NUnequalToZero[x_] := Boolean[BooleNUnequalToZero[x]]
				NNegative[x_] := Boolean[BooleNNegative[x]]
				NNonNegative[x_] := Boolean[BooleNNonNegative[x]]
			(* N | All | Sign *)
				AllNEqualToZero[x_]  := AllBoole[BooleNEqualToZero[x]]
				AllNNegative[x_]  := AllBoole[BooleNNegative[x]]
				AllNNonNegative[x_]  := AllBoole[BooleNNonNegative[x]]
				AllNPositive[x_]  := AllBoole[BooleNPositive[x]]
				AllNNonPositive[x_]  := AllBoole[BooleNonPositive[x]]
				AllNUnequalToZero[x_]  := AllBoole[BooleNUnequalToZero[x]]
			(* N | Any | Sign *)
				AnyNEqualToZero[x_] := AnyBoole[BooleNEqualToZero[x]]
				AnyNNegative[x_] := AnyBoole[BooleNNegative[x]]
				AnyNNonNegative[x_] := AnyBoole[BooleNNonNegative[x]]
				AnyNPositive[x_] := AnyBoole[BooleNPositive[x]]
				AnyNNonPositive[x_]  := AnyBoole[BooleNNonPositive[x]]
				AnyNUnequalToZero[x_] := AnyBoole[BooleNUnequalToZero[x]]
			(* N | None | Sign *)
				NoneNEqualToZero[x_] := NoneBoole[NEqualToZero[x]]
				NoneNNegative[x_] := NoneBoole[NNegative[x]]
				NoneNNonNegative[x_] := NoneBoole[NNonNegative[x]]
				NoneNPositive[x_] := NoneBoole[NPositive[x]]
				NoneNNonPositive[x_] := NoneBoole[NNonPositive[x]]
				NoneNUnequalToZero[x_] := NoneBoole[NUnequalToZero[x]]
		(* N | Compare *)
			NEqual[x_,y_] := AllNEqual[x,y]
			NGreater[x_,y_] := AllNGreater[x,y]
			NGreaterEqual[x_,y_] := AllNGreaterEqual[x,y]
			NLess[x_,y_] := AllNLess[x,y]
			NLessEqual[x_,y_] := AllNLessEqual[x,y]
			NUnequal[x_,y_] := AllNUnequal[x,y]
			(* N | Compare | Operators *)
				NEqualTo[y_][x_] := NEqual[x,y]
				NGreaterThan[y_][x_] := NGreater[x,y]
				NGreaterEqualThan[y_][x_] := NGreaterEqual[x,y]
				NLessThan[y_][x_] := NLess[x,y]
				NLessEqualThan[y_][x_] := NLessEqual[x,y]
				NUnequalTo[y_][x_] := NUnequal[x,y]
			(* N | List | Compare *)
				ListNEqual[x_,y_] := Boolean[BooleNEqual[x,y]]
				ListNGreater[x_,y_] := Boolean[BooleNGreater[x,y]]
				ListNGreaterEqual[x_,y_] := Boolean[BooleNGreaterEqual[x,y]]
				ListNLess[x_,y_] := Boolean[BooleNLess[x,y]]
				ListNLessEqual[x_,y_] := Boolean[BooleNLessEqual[x,y]]
				ListNUnequal[x_,y_] := Boolean[BooleNUnequal[x,y]]
				(* operators *)
					ListNEqualTo[y_][x_] := ListNEqual[x,y]
					ListNGreaterThan[y_][x_] := ListNGreater[x,y]
					ListNGreaterEqualThan[y_][x_] := ListNGreaterEqual[x,y]
					ListNLessThan[y_][x_] := ListNLess[x,y]
					ListNLessEqualThan[y_][x_] := ListNLessEqual[x,y]
					ListNUnequalTo[y_][x_] := ListNUnequal[x,y]
			(* N | All | Compare *)
				AllNEqual[x_,y_] := AllBoole[BooleNEqual[x,y]]
				AllNGreater[x_,y_] := AllBoole[BooleNGreater[x,y]]
				AllNGreaterEqual[x_,y_] := AllBoole[BooleNGreaterEqual[x,y]]
				AllNLess[x_,y_] := AllBoole[BooleNLess[x,y]]
				AllNLessEqual[x_,y_] := AllBoole[BooleNLessEqual[x,y]]
				AllNUnequal[x_,y_] := AllBoole[BooleNUnequal[x,y]]
				(* operators *)
					AllNEqualTo[y_][x_] := AllNEqual[x,y]
					AllNGreaterThan[y_][x_] := AllNGreater[x,y]
					AllNGreaterEqualThan[y_][x_] := AllNGreaterEqual[x,y]
					AllNLessThan[y_][x_] := AllNLess[x,y]
					AllNLessEqualThan[y_][x_] := AllNLessEqual[x,y]
					AllNUnequalTo[y_][x_] := AllNUnequal[x,y]
			(* N | Any | Compare *)
				AnyNEqual[x_,y_] := AnyBoole[BooleNEqual[x,y]]
				AnyNGreater[x_,y_] := AnyBoole[BooleNGreater[x,y]]
				AnyNGreaterEqual[x_,y_] := AnyBoole[BooleNGreaterEqual[x,y]]
				AnyNLess[x_,y_] := AnyBoole[BooleNLess[x,y]]
				AnyNLessEqual[x_,y_] := AnyBoole[BooleNLessEqual[x,y]]
				AnyNUnequal[x_,y_] := AnyBoole[BooleNUnequal[x,y]]
				(* operators *)
					AnyNEqualTo[y_][x_] := AnyNEqual[x,y]
					AnyNGreaterThan[y_][x_] := AnyNGreater[x,y]
					AnyNGreaterEqualThan[y_][x_] := AnyNGreaterEqual[x,y]
					AnyNLessThan[y_][x_] := AnyNLess[x,y]
					AnyNLessEqualThan[y_][x_] := AnyNLessEqual[x,y]
					AnyNUnequalTo[y_][x_] := AnyNUnequal[x,y]
			(* N | None | Compare *)
				NoneNEqual[x_,y_] := NoneBoole[BooleNEqual[x,y]]
				NoneNGreater[x_,y_] := NoneBoole[BooleNGreater[x,y]]
				NoneNGreaterEqual[x_,y_] := NoneBoole[BooleNGreaterEqual[x,y]]
				NoneNLess[x_,y_] := NoneBoole[BooleNLess[x,y]]
				NoneNLessEqual[x_,y_] := NoneBoole[BooleNLessEqual[x,y]]
				NoneNUnequal[x_,y_] := NoneBoole[BooleNUnequal[x,y]]
				(* operators *)
					NoneNEqualTo[y_][x_] := NoneNEqual[x,y]
					NoneNGreaterThan[y_][x_] := NoneNGreater[x,y]
					NoneNGreaterEqualThan[y_][x_] := NoneNGreaterEqual[x,y]
					NoneNLessThan[y_][x_] := NoneNLess[x,y]
					NoneNLessEqualThan[y_][x_] := NoneNLessEqual[x,y]
					NoneNUnequalTo[y_][x_] := NoneNUnequal[x,y]


(* ::Subsection::Closed:: *)
(*Vector*)


	(* Vector *)
		(* Array *)
			(* Array | Transpose *)
				ArrayTranspose[array_?ArrayQ] := Transpose[array,TwoWayRule[1,ArrayDepth[array]]]

			(* Array | Dot *)
				ArrayDot[___,{},___] := {}
				ArrayDot[a_?ArrayQ,b_?ArrayQ] := Dot[a,b]
				ArrayDot[A_?ListQ,b_?ArrayQ] := Map[Function[ArrayDot[#,b]],A]
				ArrayDot[a_?ArrayQ,B_?ListQ] := Map[Function[ArrayDot[a,#]],B]
				ArrayDot[a_,b_,c__] := ArrayDot[ArrayDot[a,b],c]

		(* Vector *)
			(* Vector | Construct *)
				OriginVector[n_Integer] := ConstantArray[0,n]
				DiagonalVector[n_] := ConstantArray[1,n]

			(* Vector | Query *)
				OriginVectorQ[vectors_] := ListNEqual[VectorNorm[vectors],0]
				UnitVectorQ[vectors_] := ListNEqual[VectorNorm[vectors],1]

		(* Dimension *)
			(* Dimension| Embedding *)
				EmbeddingDimension[vertices_] := If[ArrayQ[vertices],Last[Dimensions[vertices]],EmbeddingDimension[First[vertices]]]
					(* alias *)
					VectorDimension[vertices_] := EmbeddingDimension[vertices]

			(* Dimension | VectorSpan *)
				VectorSpanDimension[{}] := -1
				VectorSpanDimension[vectors_List] := MatrixRank[N[vectors]]
					(* alias *)
					VectorSpaceDimension[vectors_] := VectorSpanDimension[vectors]

			(* Dimension | AffineSpan *)
				AffineSpanDimension[ {} | {{}} ] := -1
				AffineSpanDimension[points_] := VectorSpanDimension[VectorTranslate[points,N[Minus[First[points]]]]]
					(* alias *)
					AffineSpaceDimension[points_] := AffineSpanDimension[points]

			(* Vector | Properties *)
				VectorCoordinateSum[vectors_] := ArrayDot[vectors,DiagonalVector[EmbeddingDimension[vectors]]]
				VectorArg[points_] := Arg[ArrayDot[VectorProject[points,2],{1,I}]]

			(* Vector | Norm *)
				VectorNormSquared[vectors_] := VectorCoordinateSum[Power[Abs[vectors],2]]
				VectorNorm[vectors_] := Sqrt[VectorNormSquared[vectors]]
				VectorNormalize[vectors_List] := Times[vectors,NInvertUnits[VectorNorm[vectors]]]

		(* Matrix *)
			(* Matrix | Query *)
				VectorIndependentQ[vectors_] := Equal[VectorSpanDimension[vectors],Length[vectors]]
				MatrixRankFullQ[vectors_] := VectorIndependentQ[vectors]

				(* alias *)
				VectorDepenentQ[vectors_] := Not[VectorIndependentQ[vectors]]
				MatrixRankMaximalQ[vectors_] := Equal[VectorSpanDimension[vectors],EmbeddingDimension[vectors]]
				MatrixRankNullQ[vectors_] := Equal[0,VectorSpanDimension[vectors]]
				MatrixInvertibleQ[matrix_] := And[SquareMatrixQ[matrix],NUnequalToZero[Det[N[matrix]]]]
					MatrixGeneralLinearQ[matrix_] := MatrixInvertibleQ[matrix]
				MatrixSpecialLinearQ[matrix_] := And[MatrixGeneralLinearQ[matrix],NEqual[Det[matrix],1]]
				NOrthogonalize[matrix_] := If[OrthogonalMatrixQ[N[matrix]], matrix, NChop[Orthogonalize[matrix]]]

			(* Matrix | Standardize *)
				MatrixSpecialize[matrix_] := MatrixNormalizeDeterminant[VectorOrient[matrix]]
				MatrixNormalizeDeterminant[matrix_] := If[MatrixGeneralLinearQ[matrix], Times[matrix,Power[Abs[Det[matrix]],Minus[Length[matrix]]]], matrix]

			(* Matrix | Random | Orthogonal *)
				RandomOrthogonalMatrix[embeddingDimension_Integer] := Orthogonalize[RandomReal[{-1,1},{embeddingDimension,embeddingDimension}]]
				RandomOrthogonalMatrix[embeddingDimension_Integer, subspaceDimension_Integer] := Orthogonalize[RandomReal[{-1,1},{subspaceDimension,embeddingDimension}]]
				RandomOrthogonalMatrix[{subspaceDimension_Integer, embeddingDimension_Integer}] := Orthogonalize[RandomReal[{-1,1},{subspaceDimension,embeddingDimension}]]

			(* Matrix | Random | SpecialOrthogonal *)
				RandomSpecialOrthogonalMatrix[dimension_Integer] := MatrixSpecialize[RandomOrthogonalMatrix[dimension]]

			(* Matrix | Random | Reflection *)
				RandomReflectionMatrix[dimension_] := ReflectionMatrix[RandomReal[{-1,1},{dimension}]]

		(* Vector | Space *)
			(* Vector | SpanningSet | Indices *)
				VectorSpanningSetIndices[{}] := {}
				VectorSpanningSetIndices[vectors_?MatrixQ, maxDimension_:Infinity] :=
					Module[{
							vectorsN,
							embeddingDimension,
							remainingVectorBooles,
							remainingVectors,
							remainingVectorIndices,
							spanningVectorIndices,
							spanningVectorIndex,
							spanningVector
						},
						embeddingDimension = EmbeddingDimension[vectors];
						vectorsN = N[vectors];
						remainingVectors = vectorsN;
						remainingVectorIndices = ListIndices[vectorsN];
						spanningVectorIndices = {};
						Do[
							remainingVectorBooles = BooleNPositive[VectorNorm[remainingVectors]];
							If[NoneBoole[remainingVectorBooles],Break[];];
							remainingVectorIndices = Pick[remainingVectorIndices, remainingVectorBooles, 1];
							remainingVectors = Pick[remainingVectors, remainingVectorBooles, 1];

							spanningVectorIndex = First[remainingVectorIndices];
							spanningVector = First[remainingVectors];
							AppendTo[spanningVectorIndices,spanningVectorIndex];
							If[EmptySetQ[remainingVectorIndices],Break[];];
							remainingVectorIndices = Rest[remainingVectorIndices];
							remainingVectors = Rest[remainingVectors];
							remainingVectors = VectorProject[remainingVectors,NullSpace[{spanningVector}]] // Chop;

						,{Min[Length[vectors],maxDimension,embeddingDimension]}];
						Return[spanningVectorIndices];
					]
			(* Vector | SpanningSet *)
				VectorSpanningSet[vectors_, maxDimension_:Infinity] := Part[vectors,VectorSpanningSetIndices[vectors,maxDimension]]
		(* Basis *)
			(* Basis | Orientation *)
				VectorOrientedQ[basis_] := If[SquareMatrixQ[basis],Positive[Chop[Det[N[basis]]]],True]
				VectorOrient[basis_] := If[Not[VectorOrientedQ[basis]],ReverseBasisOrientation[basis],basis]
			(* Basis | Orientation | Reverse *)
				ReverseBasisOrientation[ {} | {{}} ] := {}
				ReverseBasisOrientation[{v1_}] := {Minus[v1]}
				ReverseBasisOrientation[{v1_,v2_,v3___}] := {v2,v1,v3}
		(* Vector | Span *)
			(* Empty | Span *)
				EmptySpace[] := EmptySpace[<|"EmbeddingDimension"->Minus[1]|>]
				EmptySpace[n_Integer] := EmptySpace[<|"EmbeddingDimension"->n|>]
				EmptySpace[data_Association]["Dimension"] := Minus[1]
				EmptySpace[data_Association]["EmbeddingDimension"] := data["EmbeddingDimension"]
			(* Vector | Space | Construct *)
				VectorSpan[n_Integer] := VectorSpan[IdentityMatrix[n]]
				VectorSpan[n_Integer, k_Integer] := VectorSpan[IdentityMatrix[{k,n}]]
				VectorSpan[{} | {{}}] := EmptySpace[]
				VectorSpan[vector_?VectorQ] := VectorSpan[{vector}]
				VectorSpan[vectors_?MatrixQ] :=
					Module[{
							embeddingDimension,
							dimension,
							spanningVectors,
							basis
						},
						embeddingDimension = EmbeddingDimension[vectors];
						spanningVectors = VectorSpanningSet[vectors];
						dimension = Length[spanningVectors];
						basis = NOrthogonalize[spanningVectors];
						nullSpaceBasis = If[Positive[dimension], NOrthogonalize[NullSpace[basis]], IdentityMatrix[embeddingDimension]];
						VectorSpan[
							<|
								"Dimension" -> dimension,
								"EmbeddingDimension" -> embeddingDimension,
								"Basis" -> basis,
								"NullSpace" -> nullSpaceBasis
							|>
						]
					]
			(* Vector | Span | Attributes *)
				VectorSpan[data_Association]["Dimension"|"dimension"|"Dim"|"dim"|"D"|"d"] := data["Dimension"]
				VectorSpan[data_Association]["EmbeddingDimension"|"embdim"|"edim"|"EDim"|"ed"|"e"] := data["EmbeddingDimension"]
				VectorSpan[data_Association]["Basis"|"basis"|"B"|"b"] := data["Basis"]
				VectorSpan[data_Association]["NullSpace"|"nullspace"|"NullBasis"|"nullbasis"|"N"|"n"] := data["NullSpace"]
			(* Vector | Span | Properties *)
				VectorSpan[data_Association]["Codimension"|"Codim"|"codim"|"cdim"|"cd"|"c"] := Subtract[data["EmbeddingDimension"],data["Dimension"]]
				VectorSpan[data_Association]["Origin"|"origin"|"O"|"o"] := OriginVector[data["EmbeddingDimension"]]
			(* Vector | Intersection *)
				VectorSpanIntersection[___, emptySpace_EmptySpace, ___] := emptySpace
				VectorSpanIntersection[a_, b_, c__] := VectorSpanIntersection[VectorSpanIntersection[a,b],c]
				VectorSpanIntersection[a_AffineSpan, b_AffineSpan] := VectorSpanIntersection[VectorSpan[a],VectorSpan[b]]
				VectorSpanIntersection[a_AffineSpan, b_VectorSpan] := VectorSpanIntersection[VectorSpan[a],b]
				VectorSpanIntersection[a_VectorSpan, b_AffineSpan] := VectorSpanIntersection[a,VectorSpan[b]]
				VectorSpanIntersection[a_VectorSpan, b_VectorSpan] :=
					Module[{
							embeddingDimension,
							dimension,
							basis1,
							basis2,
							basis,
							nullSpaceBasis
						},
						If[Unequal[a["EmbeddingDimension"],b["EmbeddingDimension"]],Return[EmptySpace[]];];
						If[Or[a["Dimension"]==0,b["Dimension"]==0],Return[VectorSpan[embeddingDimension,0]];];
						embeddingDimension = a["EmbeddingDimension"];

						If[And[a["Dimension"]==embeddingDimension,b["Dimension"]==embeddingDimension],Return[VectorSpan[embeddingDimension,embeddingDimension]];];
						If[And[a["Dimension"]==embeddingDimension],Return[b];];
						If[And[b["Dimension"]==embeddingDimension],Return[a];];
						basis = NOrthogonalize[NullSpace[Join[a["NullSpace"],b["NullSpace"]]]];
						nullSpaceBasis = If[basis=={},
							IdentityMatrix[embeddingDimension],
							Orthogonalize[NullSpace[basis]]
						];
						dimension = Length[basis];
						VectorSpan[
							<|
								"Dimension" -> dimension,
								"EmbeddingDimension" -> embeddingDimension,
								"Basis" -> basis,
								"NullSpace" -> nullSpaceBasis
							|>
						]
					]
		(* Affine *)
			(* Affine | SpanningSet *)
				AffineSpanningSet[ {} | {{}} ] := {}
				AffineSpanningSet[point_?VectorQ] := {point}
				AffineSpanningSet[{point_?VectorQ}] := {point}
				AffineSpanningSet[points_?MatrixQ] := Prepend[Part[{points},VectorSpanningSetIndices[VectorTranslate[{points},Minus[point]]]], point]

			(* Affine | Span *)
				AffineSpan[n_Integer] := AffineSpan[n,n]
				AffineSpan[n_Integer, k_Integer] := AffineSpan[OriginVector[n],IdentityMatrix[{k,n}]]
				AffineSpan[{} | {{}}] := EmptySpace[]
				AffineSpan[point_?VectorQ, k_Integer] := AffineSpan[point,IdentityMatrix[{k,EmbeddingDimension[point]}]]
				AffineSpan[point_?VectorQ] := AffineSpan[point,{}]
				AffineSpan[{point_?VectorQ}] := AffineSpan[point,{}]
				AffineSpan[{point_?VectorQ, points__?VectorQ}] := AffineSpan[point, VectorTranslate[{points},Minus[point]]]

				AffineSpan[point_?VectorQ, {} | {{}} ] := AffineSpan[point,0]
				AffineSpan[point_?VectorQ, vector_?VectorQ] := AffineSpan[point, {vector}]
				AffineSpan[point_?VectorQ, vectors_?MatrixQ] :=
					Module[{
							embeddingDimension,
							spanningVectors,
							affineDimension,
							affineBasis,
							nullSpaceBasis,
							affineOrigin
						},
						embeddingDimension = EmbeddingDimension[point];
						spanningVectors = VectorSpanningSet[vectors];
						affineDimension = Length[spanningVectors];
						affineOrigin = point;
						affineBasis = spanningVectors // NOrthogonalize;
						nullSpaceBasis = NullSpace[Prepend[affineBasis,OriginVector[embeddingDimension]]] // NOrthogonalize;
						AffineSpan[
							<|
								"Dimension" -> affineDimension,
								"EmbeddingDimension" -> embeddingDimension,
								"Basis" -> affineBasis,
								"NullSpace" -> nullSpaceBasis,
								"Origin" -> affineOrigin
							|>
						]
					]
			(* AffineSpan | Attributes *)
				AffineSpan[data_Association]["Dimension"|"dimension"|"Dim"|"dim"|"D"|"d"] := data["Dimension"]
				AffineSpan[data_Association]["EmbeddingDimension"|"embdim"|"edim"|"EDim"|"ed"|"e"] := data["EmbeddingDimension"]
				AffineSpan[data_Association]["Basis"|"basis"|"B"|"b"] := data["Basis"]
				AffineSpan[data_Association]["NullSpace"|"nullspace"|"NullBasis"|"nullbasis"|"N"|"n"] := data["NullSpace"]
				AffineSpan[data_Association]["Origin"|"origin"|"O"|"o"] := data["Origin"]
			(* AffineSpan | Properties *)
				AffineSpan[data_Association]["Codimension"|"Codim"|"codim"|"cdim"|"cd"|"c"] := Subtract[data["EmbeddingDimension"],data["Dimension"]]
			(* ->> EmptySpace *)
				EmptySpace[AffineSpan[data_Association]] := EmptySpace[data["EmbeddingDimension"]]
				EmptySpace[VectorSpan[data_Association]] := EmptySpace[data["EmbeddingDimension"]]
				VectorSpan[EmptySpace[data_Association]] := EmptySpace[data["EmbeddingDimension"]]
				AffineSpan[EmptySpace[data_Association]] := EmptySpace[data["EmbeddingDimension"]]
			(* ->> AffineSpan *)
				AffineSpan[VectorSpan[data_Association]] := AffineSpan[Append[data,"Origin"->OriginVector[data["EmbeddingDimension"]]]]
			(* ->> VectorSpan *)
				VectorSpan[AffineSpan[data_Association]] := VectorSpan[Delete[data,"Origin"]]
			(* ->> AffineSpace ->> Graphics(3D) *)
				AffineSpan /: AffineSpace[AffineSpan[data_Association]] := AffineSpace[data["Origin"],data["Basis"]]
				VectorSpan /: AffineSpace[VectorSpan[data_Association]] := AffineSpace[OriginVector[data["EmbeddingDimension"]],data["Basis"]]
				EmptySpace /: AffineSpace[EmptySpace[data_Association]] := EmptyRegion[data["EmbeddingDimension"]]
			(* AffineSpanIntersection *)
				AffineSpanIntersection[a_,b_,c__] := AffineSpanIntersection[AffineSpanIntersection[a,b],c]
				AffineSpanIntersection[___, EmptySpace[n_], ___] := EmptySpace[n]
				AffineSpanIntersection[vec1_VectorSpan, vec2_VectorSpan] := AffineSpanIntersection[AffineSpan[vec1], AffineSpan[vec2]]
				AffineSpanIntersection[vec_VectorSpan, aff_AffineSpan] := AffineSpanIntersection[aff, AffineSpan[vec]]
				AffineSpanIntersection[aff_AffineSpan, vec_VectorSpan] := AffineSpanIntersection[aff, AffineSpan[vec]]
				AffineSpanIntersection[aff1_AffineSpan, aff2_AffineSpan] :=
					Module[{
							embeddingDimension,
							dimension1,
							dimension2,
							point1,
							basis1,
							point2,
							basis2,
							intersectionBasis,
							liftedBasis1,
							liftedBasis2,
							liftedIntersectionBasis,
							liftedIntersectionPoint,
							intersectionPoint
						},
						If[Unequal[aff1["EmbeddingDimension"],aff2["EmbeddingDimension"]],
							Message[VectorSpanIntersection::DimensionMismatch];
							Return[AffineSpan[]];
						];
						embeddingDimension = aff1["EmbeddingDimension"];
						dimension1 = aff1["Dimension"];
						dimension2 = aff2["Dimension"];
						point1 = aff1["Origin"];
						basis1 = aff1["Basis"];
						point2 = aff2["Origin"];
						basis2 = aff2["Basis"];
						zeroBasis = ConstantArray[0,{1,embeddingDimension}];
						If[basis1 == {}, basis1 = zeroBasis;];
						If[basis2 == {}, basis2 = zeroBasis;];
						intersectionBasis = NOrthogonalize[NullSpace[Join[NullSpace[basis1],NullSpace[basis2],zeroBasis]]];
						lisftedZeroBasis = ConstantArray[0,{1,Plus[embeddingDimension,1]}];
						liftedBasis1 = ArrayFlatten[{{basis1,0},{{point1},1}}];
						liftedBasis2 = ArrayFlatten[{{basis2,0},{{point2},1}}];
						liftedIntersectionBasis = NullSpace[Join[NullSpace[liftedBasis1],NullSpace[liftedBasis2],lisftedZeroBasis]];
						If[Or[liftedIntersectionBasis=={},AllNEqualToZero[Part[liftedIntersectionBasis,All,-1]]],
							Return[EmptyRegion[embeddingDimension]];];
						liftedIntersectionPoint = SelectFirst[liftedIntersectionBasis, Last/*Chop/*UnequalTo[0]];
						liftedIntersectionPoint = Divide[liftedIntersectionPoint,Last[liftedIntersectionPoint]];
						intersectionPoint = Most[liftedIntersectionPoint];
						intersectionPoint = intersectionBasis . intersectionPoint;
						Return[AffineSpan[intersectionPoint,intersectionBasis]];
					]

		(* Vector *)
			(* Vector | Translate *)
				VectorTranslate[{},vector_] := {}
				VectorTranslate[points_,{}] := {}
				VectorTranslate[points_?ArrayQ,vector_?VectorQ] := ArrayTranspose[Plus[ArrayTranspose[points],vector]]
				VectorTranslate[points_?ListQ,vector_?VectorQ] := Map[Function[VectorTranslate[#,vector]],points]
				VectorTranslate[points_?ListQ,vectors_?ListQ] := Map[Function[VectorTranslate[points,#]],vectors]
				VectorTranslate[points_,value_?(Composition[Not,ListQ])] := Plus[points,value]
				VectorTranslate[vectors_][points_] := VectorTranslate[points,vectors]

			(* Vector | Centroid *)
				VectorCentroid[point_?VectorQ] := point
				VectorCentroid[points_?MatrixQ] := Mean[points]
				VectorCentroid[points_List] := Map[VectorCentroid,points]

			(* Vector | Recenter *)
				VectorRecenter[point_?VectorQ] := point
				VectorRecenter[points_?MatrixQ] := VectorTranslate[points,Minus[Mean[points]]]
				VectorRecenter[points_List] := Map[VectorRecenter,points]

			(* Vector | Dilate *)
				VectorDilate[points_, scale_, pivot_] := VectorTranslate[Times[VectorTranslate[points,Minus[pivot]],scale],Plus[pivot]]
				VectorDilate[points_?MatrixQ, scale_] := VectorDilate[points,scale,Mean[points]]
				VectorDilate[points_List, scale_] := Map[Function[VectorDilate[#,scale]],points]
				VectorDilate[scale_, pivot_?VectorQ][points_] := VectorDilate[points,scale,pivot]
				VectorDilate[scale_][points_] := VectorDilate[points,scale]

			(* Vector | Sort *)
				VectorSort[vector_?VectorQ] := vector
				VectorSort[vectors_?MatrixQ] := NumericalSort[vectors]
				VectorSort[vectors_List] := NumericalSort[Map[VectorSort,vectors]]

			(* Vector | BoundedQ *)
				VectorBoundedQ[points_, normals_, bounds_] := Boolean[VectorBoundedQBoole[points,normals,bounds]]
				VectorBoundedQ[normals_, bounds_][points_] := VectorBoundedQ[points,normals,bounds]

			(* Vector | BoundedQ | Boole*)
				VectorBoundedQBoole[points_, normal_?VectorQ, bound_] := BooleNNonPositive[Subtract[ArrayDot[points,normal],bound]]
				VectorBoundedQBoole[points_, normals_?MatrixQ, bounds_] := 	BooleEqual[ArrayDot[BooleNNonPositive[VectorTranslate[ArrayDot[points,Transpose[normals]],-bounds]],DiagonalVector[Length[normals]]],Length[normals]]
				VectorBoundedQBoole[normals_, bounds_][points_] := VectorBoundedQBoole[points,normals,bounds]

			(* Vector | Match *)
				VectorMatchQ[pointsA_,pointsB_] := NEqual[VectorSort[pointsA],VectorSort[pointsB]]

			(* Vector | Cases *)
				VectorCases[points_?MatrixQ, pattern_?MatrixQ] :=
					Module[{
							patternPoints,
							patternVectors,
							patternVectorIndex,
							patternOffsetVector,
							patternVector,
							matchedPatternBasePoints
						},
						patternPoints = DeleteDuplicates[pattern];
						patternBasePoint = First[patternPoints];
						patternVectors = VectorTranslate[Rest[patternPoints],Minus[patternBasePoint]];
						matchedPatternBasePoints = points;
						Monitor[
							Do[
								patternVector = Part[patternVectors,patternVectorIndex];
								matchedPatternBasePoints = Intersection[matchedPatternBasePoints,VectorTranslate[points,Minus[patternVector]]];
							,{patternVectorIndex,Length[patternVectors]}];
						,ProgressIndicator[patternVectorIndex,{1,Length[patternVectors]}]
						,1
						];
						matchedPatternBasePoints = VectorTranslate[matchedPatternBasePoints,Minus[patternBasePoint]];
						Return[matchedPatternBasePoints];
					]

			(* N | Vector | Cases *)
				NVectorCases[points_?MatrixQ, pattern_?MatrixQ] :=
					Module[{
							pointsN,
							patternPoints,
							patternVectors,
							patternVectorIndex,
							patternOffsetVector,
							patternVector,
							matchedPatternBasePointsN,
							matchedPatternBasePoints
						},
						pointsN = N[points];
						patternPoints = NDeleteDuplicates[pattern];
						patternBasePoint = First[patternPoints];
						patternVectors = VectorTranslate[Rest[patternPoints],Minus[patternBasePoint]];
						matchedPatternBasePointsN = NRound[pointsN];
						Monitor[
							Do[
								patternVector = Part[patternVectors,patternVectorIndex];
								matchedPatternBasePointsN = Intersection[matchedPatternBasePointsN,NRound[VectorTranslate[pointsN,Minus[patternVector]]]];
							,{patternVectorIndex,Length[patternVectors]}];
						,ProgressIndicator[patternVectorIndex,{1,Length[patternVectors]}]
						,1
						];
						matchedPatternBasePoints = NIntersection[points,matchedPatternBasePointsN];
						matchedPatternBasePoints = VectorTranslate[matchedPatternBasePoints,Minus[patternBasePoint]];
						Return[matchedPatternBasePoints];
					]
			(* Vector | Shells *)
				VectorShells[points_] := SortBy[Parts[points,GatherIndices[NRound[VectorNorm[N[points]]]]],Composition[Norm,N,First]]

			(* Vector | Slice | Points *)
				VectorSlicePoints[{},normal_,value_] := {}
				VectorSlicePoints[points_,normal_,value_] := Pick[points, BooleNEqual[ArrayDot[N[points],N[normal]],N[value]], 1]
				VectorSlicePoints[normal_,value_][points_] := VectorSlicePoints[points,normal,value]

			(* Vector | Slice | Points | Maximal *)
				VectorMaximalByNormal[points_,normal_] := Part[points,NMaximalIndices[Dot[N[points],normal]]]

			(* Vector | Chop | Points *)
				VectorChopPoints[{},normal_,value_] := {}
				VectorChopPoints[points_?ArrayQ,normal_,value_] := Pick[points,NonPositive[Sign[Chop[Subtract[Dot[points,normal],value]]]]]
				VectorChopPoints[points_,normal_,value_] := Map[Function[VectorChopPoints[#,normal,value]],points]
				VectorChopPoints[normal_,value_][points_] := VectorChopPoints[points,normal,value]

			(* Vector | Transform *)
				VectorTransform[points_,matrix_] := ArrayDot[points,Transpose[matrix]]
				VectorTransform[matrix_][points_] := VectorTransform[points,matrix]

			(* Vector | Rotate *)
				VectorRotate[points_,matrix_?MatrixQ] := VectorTransform[points,matrix]
				VectorRotate[points_,theta_?NumericQ] := VectorTransform[points,RotationMatrix[theta,IdentityMatrix[{2,EmbeddingDimension[points]}]]]
				VectorRotate[matrix_][points_] := VectorTransform[points,matrix]

			(* Vector | Affine | Rotate *)
				VectorAffineRotate[points_,transform_,pivot_] := VectorTranslate[VectorRotate[VectorTranslate[points,Minus[pivot]],transform],pivot]
				VectorAffineRotate[transform_,pivot_][points_] := VectorAffineRotate[points,transform,pivot]

			(* Vector | Projection *)
				VectorProject[{},subspace_] := {}
				VectorProject[points_,0] := ArrayDot[points,ConstantArray[{},EmbeddingDimension[points]]]
				VectorProject[points_,dimension_Integer] := ArrayDot[points,IdentityMatrix[{EmbeddingDimension[points],dimension}]]
				VectorProject[points_,{}] := ArrayDot[points,ConstantArray[{},EmbeddingDimension[points]]]
				VectorProject[points_,vector_?VectorQ] := ArrayDot[points,vector]
				VectorProject[points_,subspace_?MatrixQ] := ArrayDot[points,Transpose[subspace]]
				VectorProject[subspace_][points_] := VectorProject[points,subspace]

			(* Vector | Affine | Projection *)
				VectorAffineProject[points_,subspace_,pivot_] := VectorProject[VectorTranslate[points,Minus[pivot]],subspace]
				VectorAffineProject[subspace_,pivot_][points_] := VectorAffineProject[points,subspace,pivot]

			(* Vector | Embed *)
				VectorEmbed[{},subspace_] := {}
				VectorEmbed[points_,dimension_Integer] := ArrayDot[points,IdentityMatrix[{EmbeddingDimension[points],dimension}]]
				VectorEmbed[points_,basis_?MatrixQ] := ArrayDot[points,basis]
				VectorEmbed[basis_][points_] := VectorEmbed[points,basis]

			(* Vector | Affine *)
				VectorAffineEmbed[points_,subspace_,pivot_] := VectorTranslate[VectorEmbed[points,subspace],pivot]
				VectorAffineTransform[points_,xform_?SquareMatrixQ,pivot_] := VectorTranslate[VectorTransform[VectorTranslate[points,Minus[pivot]],xform],pivot]

			(* Vector | MatrixPlot *)
				VectorMatrixPlot[vectors_?MatrixQ] := MatrixPlot[Transpose[vectors], Rule[Frame,None], Rule[Mesh,All], Rule[MaxPlotPoints,Infinity], Rule[MeshStyle,{{AbsoluteThickness[1],Hue[0,0,0,.1]}}], Rule[PixelConstrained,{8,8}]]
				VectorMatrixPlot[vector_?VectorQ] := VectorMatrixPlot[{vector}]
				VectorMatrixPlot[matrices_?ListQ] := Map[VectorMatrixPlot,matrices]

			(* Vector Distance *)
				(* Vector | Distance | point *)
					VectorDistance[points_,origin_?Vector] := VectorNorm[VectorTranslate[points,Minus[origin]]]
				(* Vector | Distance | Line *)
					VectorDistanceToLine[base_,axis_][points_] := VectorNorm[VectorProject[VectorTranslate[points,Minus[base]],Orthogonalize[NullSpace[{axis}]]]]

			(* Path *)
				(* Vector | Path | Edges *)
					VectorPathEdges[path_?MatrixQ] := Partition[path,2,1]

				(* Vector | Path | ParallelTransport *)
					VectorPathParallelTransport[path_] :=
						Module[{
								dim,
								rotations,
								frames,
								tangent,
								tangentNew,
								frame,
								rotation
							},
							dim = EmbeddingDimension[path];
							frame = IdentityMatrix[dim];
							rotations = {};
							frame = IdentityMatrix[dim];
							frames = {};
							tangent = First[frame];
							Do[
								tangentNew = VectorNormalize[path[[i+1]]-path[[i]]];
								If[MatrixRank[{tangent,tangentNew}] == 2,
									rotation = RotationMatrix[{tangent,tangentNew}];
									,
									If[Negative[Dot[tangent,tangentNew]],
										rotation = ReflectionMatrix[tangentNew];
										,
										rotation = IdentityMatrix[dim];
									];
								];
								tangent = tangentNew;
								frame = Dot[frame,Transpose[rotation]];
								AppendTo[frames,frame];
								AppendTo[rotations,rotation];
							,{i,Length[path]-1}
							];
							AppendTo[rotations,Last[rotations]];
							AppendTo[frames,Last[frames]];
							Return[frames];
						]

			(* Simplex *)
				(* Simplex | Normal *)
					SimplexNormal[{{x_}}] := {1}
					SimplexNormal[simplex_] := Apply[Cross,VectorTranslate[Rest[simplex],Minus[First[simplex]]]]

				(* Simplex | Measure *)
					SimplexMeasure /: SimplexMeasure[Simplex[pointCoordinates_]] := Abs[SimplexMeasure[pointCoordinates]]
					SimplexMeasure[simplexInput_?MatrixQ] :=
						Module[{
								simplexVertices,
								embeddingDimension,
								affineDimension,
								simplexVectors,
								measure
							},
							simplexVertices = simplexInput;
							embeddingDimension = EmbeddingDimension[simplexVertices];
							affineDimension = Subtract[Length[simplexVertices],1];
							If[Greater[affineDimension,embeddingDimension],
								Return[0];];
							simplexVectors = VectorTranslate[Rest[simplexVertices],Minus[First[simplexVertices]]];
							If[Less[affineDimension,embeddingDimension],
								simplexVectors = VectorProject[simplexVectors,Orthogonalize[simplexVectors]];];
							measure = Divide[Det[simplexVectors],Factorial[affineDimension]];
							measure = Abs[measure];
							Return[measure];
						]

				(* Simplex | OrientedFacets *)
					SimplexOrientedFacets[simplexVertices_] :=
						If[LessEqual[Length[simplexVertices],2],
							Subsets[simplexVertices,{Length[simplexVertices]-1}]
							,
							Module[{simplexVertexIndices,simplexVertexCount,simplexFacetCount,orientedFacetVertexIndices,orientedFacetVertices},
								simplexVertexIndices = ListIndices[simplexVertices];
								simplexVertexCount = Length[simplexVertices];
								simplexFacetCount = simplexVertexCount;
								orientedFacetVertexIndices = ConstantArray[Most[simplexVertexIndices],simplexFacetCount];
								orientedFacetVertexIndices[[1,{1,2}]] = orientedFacetVertexIndices[[1,{2,1}]];
								Do[
									orientedFacetVertexIndices[[index+1,index]] = simplexVertexCount;
								,{index,1,simplexVertexCount-1}];
								orientedFacetVertices = Parts[simplexVertices,orientedFacetVertexIndices];
								Return[orientedFacetVertices];
							]
						]

				(* Simplex | OrientedQ *)
					SimplexOrientedQ[simplex_?SquareMatrixQ] := Positive[Det[VectorTranslate[Rest[simplex],N[Minus[First[simplex]]]]]]
					SimplexRegularQ[ verticesInput_?MatrixQ ] :=
						Module[{
								dimension,
								vertexCount,
								edgeCount,
								vertices,
								simplexRegularQ,
								edgeLengthFunction,
								edgeLengthTarget,
								edgeIndex,
								edgeVertices,
								edgeLength
							},
							dimension = EmbeddingDimension[verticesInput];
							vertexCount = Length[verticesInput];
							edgeCount = Binomial[vertexCount,2];
							If[Equal[vertexCount,1], Return[True];];
							If[Greater[vertexCount,dimension], Return[False];];
							vertices = N[verticesInput];
							edgeLengthFunction = RightComposition[Differences,First,Norm];
							edgeLengthTarget = First[Subsets[vertices,{2},{1}]] // edgeLengthFunction;
							simplexRegularQ = True;
							Do[
								edgeVertices = First[Subsets[vertices,{2},{edgeIndex}]];
								edgeLength = edgeVertices // edgeLengthFunction;
								If[Not[NEqual[edgeLength,edgeLengthTarget]],
									simplexRegularQ = False;
									Break[];
								];
							,{edgeIndex,edgeCount}];
							Return[simplexRegularQ];
						]

			(* Polygon *)
				(* Polygon | Basis *)
					PolygonBasis[points_] := Orthogonalize[VectorSpanningSet[VectorTranslate[Rest[N[points]],Minus[First[N[points]]]],2]]

				(* Polygon | VertexOrdering *)
					PolygonVertexOrdering[pointsInput_] :=
						Module[{
								points,
								vectors,
								plane,
								args,
								norms,
								ordering
							},
							points = N[pointsInput];
							If[LessEqual[Length[points],3], Return[ListIndices[points]];];
							vectors = VectorRecenter[points];
							plane = Orthogonalize[VectorSpanningSet[vectors,2]];
							points = VectorProject[points,plane];
							ordering = Ordering[VectorArg[VectorProject[vectors,plane]]];
							Return[ordering];
						]

				(* Polygon | VertexSort *)
					PolygonVertexSort[polygon_] := Part[polygon,PolygonVertexOrdering[polygon]]

				(* Polygon | Attributes *)
					PolygonEdges[polygon_] := Transpose[{polygon,RotateLeft[polygon]}]
					PolygonMeasure[polygon_] := PolygonArea[polygon]
					PolygonArea[polygonInput_] :=
						Module[{
								polygon,
								basis,
								origin,
								polygonArea
							},
							polygon = VectorRecetner[N[polygonInput]];
							basis = Orthogonalize[VectorSpanningSet[VectorRecetner[N[polygon]],2]];
							polygon = polygon // VectorProject[basis];
							polygon = Part[polygon,VectorArg[polygon]];
							origin = OriginVector[EmbeddingDimension[polygon]];
							polygonArea = polygon // VectorPathEdges // Map[Prepend[origin]] // Map[SimplexMeasure] // Total;
							Return[polygonArea];
						]

			NOrthogonalMatrixQ[matrix_?MatrixQ] := OrthogonalMatrixQ[N[matrix]];

			PrincipalVectorAngles[basis1_?(N /* OrthogonalMatrixQ), basis2_?(N /* OrthogonalMatrixQ)] := Re@ArcCos[SingularValueList[Dot[basis1,Transpose[basis2]]]]
			PrincipalVectorAngles[basis1_?MatrixQ, basis2_?MatrixQ] := Re@ArcCos[SingularValueList[Dot[Orthogonalize[basis1],Transpose[Orthogonalize[basis2]]]]]
			PrincipalVectorAngles[a:_VectorSpan|_AffineSpan, b:_VectorSpan| q_AffineSpan] := PrincipalVectorAngles[a["Basis"],b["Basis"]]


(* ::Subsection:: *)
(*End*)


	End[];


EndPackage[];
