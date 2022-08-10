---
name: SpatioTemporal
topic: Handling and Analyzing Spatio-Temporal Data
maintainer: Edzer Pebesma, Roger Bivand
email: edzer.pebesma@uni-muenster.de
version: 2022-08-10
source: https://github.com/cran-task-views/SpatioTemporal/
---


This task view aims at presenting R packages that are useful for the
analysis of spatio-temporal data.

If something is inaccurate or missing, please send an e-mail to the
maintainer or submit an issue or pull request in the GitHub
repository linked above.

The following people contributed to this task view: Roger Bivand, Achim
Zeileis, Michael Sumner, Ping Yang.

Although one could argue that all data are spatio-temporal, as they must
have been taken somewhere and at some point in time, in many cases the
spatial locations or times of observation are not registered, and
irrelevant to the purpose of the study. Here, we will address the cases
where both location *and* time of observation are registered, and
relevant for the analysis of the data. The
`r view("Spatial")` and `r view("TimeSeries")`
task views shed light on spatial, and temporal data handling and
analysis, individually.

### Representing data

-   **In long tables:** In some cases, spatio-temporal data can be held
    in tables ( `data.frame` objects), with longitude, latitude and time
    as three of the columns, or an identifier for a location or region
    and time as columns. For instance, data sets in package
    `r pkg("plm")` for linear panel models have repeated
    observations for observational units, where these units often refer
    to spatial areas (countries, states) by an index. This index (a
    name, or number) can be matched to the spatial coordinates
    (polygons) of the corresponding area, an example of this is given by
    [Pebesma (2012, Journal of Statistical
    Software)](http://www.jstatsoft.org/v51/i07/) . As these data sets
    usually contain more than one attribute, to hold the data in a
    two-dimensional table a *long table* form is chosen, where each
    record contains the index of the observational unit, observation
    time, and all attributes.
-   **In time-wide tables:** When a single attribute is considered,
    another layout is that of the *time-wide table* , where each
    observational unit forms a record and each column an observation
    time. `r pkg("googleVis")` lets you analyze such data in
    a way similar to gapminder (see links).
-   **In space-wide tables:** An example of a space-wide table is the
    Irish wind data set, obtained by `data(wind)` in package
    `r pkg("gstat", priority = "core")`. It has time series
    as different columns, each column representing one location (weather
    station). The `stConstruct` function in package
    `r pkg("spacetime", priority = "core")` accepts data in
    long, time-wide or space-wide tables.
-   **Generic classes:** Formal classes for spatio-temporal data in R
    are provided by the `r pkg("spacetime")` package, which
    offers S4 classes for full space-time grids (every observational
    unit contains an observation for each observation time), sparse
    space-time grids (regular, but incomplete grids), irregular
    space-time data (each observational unit is observed at its own
    time), and has limited support for trajectory data.
    `r pkg("spacetime")` classes have
    `r pkg("sp", priority = "core")` and
    `r pkg("xts", priority = "core")` objects as slots for
    the spatial and temporal components, and can deal with all spatial
    classes (points, lines, polygons, grids) of
    `r pkg("sp")`, regular and irregular time series, and
    extend the powerful methods (selection, aggregation, plotting
    coercion) from both packages. Package `r pkg("stars", priority = "core")`
    replaces `r pkg("spacetime")` for full space-time grids; package
    `r github("r-spatial/sftime")` replaces its irregulat spacetime data,
    although `sf` objects from package `r pkg("sf")` can also do that.
-   **Dedicated classes:** dedicated classes are offered for:
    -   **data cubes:** package `r pkg("stars", priority = "core")` provides
        methods to create, analyse and visualise raster and vector data cubes,
        array data with arbitrary dimensionaly where some dimensions refer to
        space (raster, vector) and/or time
    -   **Gridded/raster data:** package
        `r pkg("raster", priority = "core")` and its successor 
        `r pkg("terra", priority = "core")` deal with sets
        of rasters (called bricks, or stacks), and a set may reflect a
        temporal sequence; both package also handle vector data
    -   **Lattice data:** package
        `r pkg("surveillance", priority = "core")` provides a
        class `sts`, which holds a `SpatialPolygonsDataFrame` slot for
        the areas, and numeric slots to define a regular time series (no
        time objects, such as `POSIXct`).
    -   **Point patterns:** Package `r pkg("spatstat")`
        provides a class `ppx` that deals spatial and temporal
        coordinate. None of the point pattern classes mentioned support
        spatial or explicit temporal reference systems.
    -   **Trajectory data:** see also the dedicated task view
        `r view("Tracking")`; Package
        `r pkg("adehabitatLT", priority = "core")` offers a
        class `ltraj` for trajectories, and methods for analyzing them;
        the packages `r pkg("move")` and
        `r pkg("trip", priority = "core")` both extend
        `r pkg("sp")` based classes for trajectories. A blog
        post on [tidy storm
        trajectories](http://r-spatial.org/r/2017/08/28/nest.html)
        points out how nested dataframes, along with geometry list
        columns of the `r pkg("sf")` package, can be used to
        model sets of trajectories, and visualise properties at the set
        level and at the level of individual fixes.

### Analyzing data

-   **Geostatistical data**
    -   `r pkg("gstat")` provides kriging, methods of
        moments variogram estimation and model fitting for a limited
        range of spatio-temporal models.
    -   `r pkg("IDE")` provides functionality for modelling
        spatio-temporal data using the integro-difference equation.
    -   `r pkg("RandomFields", priority = "core")` provides
        kriging, conditional simulation, and covariance functions and
        maximum likelihood function fitting for a very wide range of
        spatio-temporal covariance models.
    -   the `r pkg("spTimer")` package is able to fit,
        spatially predict and temporally forecast large amounts of
        space-time data using Bayesian Gaussian Process (GP) Models,
        Bayesian Auto-Regressive (AR) Models, and Bayesian Gaussian
        Predictive Processes (GPP) based AR Models.
    -   `r pkg("spBayes")` provides functions for fitting
        Bayesian dynamic space-time regression models for settings where
        space is viewed as continuous but time is taken to be discrete.
    -   `r pkg("Stem")` provides estimation of the
        parameters of a spatio-temporal model using the EM algorithm,
        estimation of the parameter standard errors using a
        spatio-temporal parametric bootstrap, spatial mapping.
    -   `r pkg("pastecs")` is a package for the regulation,
        decomposition and analysis of space-time series.
    -   `r pkg("STMedianPolish")` analyses spatio-temporal
        data, decomposing data in n-dimensional arrays and using the
        median polish technique.
    -   R-Forge package `r rforge("spcopula")` provides a
        framework to analyze via copulas spatial and spatio-temporal
        data provided in the format of the spacetime package.
        Additionally, support for calculating different multivariate
        return periods is implemented.
    -   `r pkg("solaR")` is a package for computing solar
        radiation and photovoltaic systems performance.
    -   `r pkg("nlme")` and `r pkg("lme4")`
        contain functions to fit linear mixed models, and have
        facilities to model spatial and/or temporal effects.
    -   `r pkg("mlr3spatiotempcv")` provides spatiotemporal resampling 
        methods; it extends the mlr3 ML framework with spatio-temporal 
        resampling methods to account for the presence of spatiotemporal 
        autocorrelation in predictor variables. 
-   **Point patterns**
    -   `r pkg("splancs")` provides methods for spatial and
        space-time point pattern analysis (khat, kernel3d, visualizing).
    -   [ptproc](http://www.biostat.jhsph.edu/~rpeng/software/)
        (off-CRAN) provides methods and classes for spatio-temporal
        ("multi-dimensional") point process.
-   **Lattice data**
    -   `r pkg("surveillance")` provides temporal and
        spatio-temporal modeling and monitoring of epidemic phenomena.
    -   `r pkg("plm")` fits linear panel models.
    -   `r pkg("splm")` provides estimation and diagnostic
        testing of econometric models for spatial panel data.
    -   `r pkg("sphet")` fit spatial models with
        heteroskedastic innovations.
    -   `r pkg("nlme")` and `r pkg("lme4")`
        contain functions to fit linear mixed models, and have
        facilities to model spatial and/or temporal effects.
    -   `r pkg("rsatscan")` provides an R interface to the
        free (but non-open source) program SaTScan.
    -   `r pkg("CARBayesST")` implements a class of
        spatio-temporal generalised linear mixed models for areal unit
        data, with inference in a Bayesian setting using Markov chain
        Monte Carlo (McMC) simulation.
    -   `r pkg("gapfill")` provides tools to fill missing
        values in satellite data and to develop new gap-fill algorithms.
        The methods are tailored to data (images) observed at
        equally-spaced points in time. The package is illustrated with
        MODIS NDVI data.
-   **Moving objects, trajectories** see the dedicated task view `r view("Tracking")`
    -   There is a large (74+) and growing number of tracking,
        trajectory, movement and related packages on CRAN. The review
        paper by [Loo et al (2018)](https://doi.org/10.1111/1365-2656.13116)
        provides a guide to summarize many available packages and their
        functionality.
    -   `r pkg("adehabitatLT")` A collection of tools for
        the analysis of animal movements.
    -   `r pkg("animalTrack")` 2D and 3D animal tracking
        data can be used to reconstruct tracks through time/space with
        correction based on known positions. 3D visualization of animal
        position and attitude.
    -   `r pkg("anipaths")` Animation of observed
        trajectories using spline-based interpolation. Intended to be
        used exploratory data analysis, and perhaps for preparation of
        presentations.
    -   `r pkg("argosfilter")` Functions to filters animal
        satellite tracking data obtained from Argos. It is especially
        indicated for telemetry studies of marine animals, where Argos
        locations are predominantly of low-quality.
    -   `r pkg("AtmRay")` Calculates acoustic traveltimes
        and ray paths in 1-D, linear atmospheres. Later versions will
        support arbitrary 1-D atmospheric models, such as radiosonde
        measurements and standard reference atmospheres.
    -   `r pkg("BayesianAnimalTracker")` Bayesian melding
        approach to combine the GPS observations and Dead-Reckoned path
        for an accurate animal's track, or equivalently, use the GPS
        observations to correct the Dead-Reckoned path. It can take the
        measurement errors in the GPS observations into account and
        provide uncertainty statement about the corrected path. The main
        calculation can be done by the BMAnimalTrack function.
    -   `r pkg("BBMM")` The model provides an empirical
        estimate of a movement path using discrete location data
        obtained at relatively short time intervals.
    -   `r pkg("bcpa")` The Behavioral Change Point Analysis
        (BCPA) is a method of identifying hidden shifts in the
        underlying parameters of a time series, developed specifically
        to be applied to animal movement data which is irregularly
        sampled. The method is based on: E. Gurarie, R. Andrews and K.
        Laidre A novel method for identifying behavioural changes in
        animal movement data (2009) Ecology Letters 12:5 395-408.
    -   `r pkg("bsam")` Tools to fit Bayesian state-space
        models to animal tracking data. Models are provided for location
        filtering, location filtering and behavioural state estimation,
        and their hierarchical versions. The models are primarily
        intended for fitting to ARGOS satellite tracking data but
        options exist to fit to other tracking data types. For Global
        Positioning System data, consider the 'moveHMM' package.
        Simplified Markov Chain Monte Carlo convergence diagnostic
        plotting is provided but users are encouraged to explore tools
        available in packages such as 'coda' and 'boa'.
    -   `r pkg("caribou")` This is a package for estimating
        the population size of migratory caribou herds based on large
        scale aggregations monitored by radio telemetry. It implements
        the methodology found in the article by Rivest et al. (1998)
        about caribou abundance estimation. It also includes a function
        based on the Lincoln-Petersen Index as applied to radio
        telemetry data by White and Garrott (1990).
    -   `r pkg("crawl")` Fit continuous-time correlated
        random walk models with time indexed covariates to animal
        telemetry data. The model is fit using the Kalman-filter on a
        state space version of the continuous-time stochastic movement
        process.
    -   `r pkg("ctmcmove")` Software to facilitates taking
        movement data in xyt format and pairing it with raster
        covariates within a continuous time Markov chain (CTMC)
        framework. As described in Hanks et al. (2015), this allows
        flexible modeling of movement in response to covariates (or
        covariate gradients) with model fitting possible within a
        Poisson GLM framework.
    -   `r pkg("ctmm")` Functions for identifying, fitting,
        and applying continuous-space, continuous-time stochastic
        movement models to animal tracking data.
    -   `r pkg("EMbC")` Unsupervised, multivariate, binary
        clustering for meaningful annotation of data, taking into
        account the uncertainty in the data. A specific constructor for
        trajectory analysis in movement ecology yields behavioural
        annotation of trajectories based on estimated local measures of
        velocity and turning angle, eventually with solar position
        covariate as a daytime indicator, ("Expectation-Maximization
        Binary Clustering for Behavioural Annotation").
    -   `r pkg("eyelinker")` Eyelink eye trackers output a
        horrible mess, typically under the form of a '.asc' file. The
        file in question is an assorted collection of messages, events
        and raw data. This R package will attempt to make sense of it.
    -   `r pkg("eyetracking")` Misc function for working
        with eyetracking data
    -   `r pkg("fishmove")` Functions to predict fish
        movement parameters plotting leptokurtic fish dispersal kernels
        (see Radinger and Wolter, 2014: Patterns and predictors of fish
        dispersal in rivers. Fish and Fisheries. 15:456-473.)
    -   `r pkg("foieGras")` Fits continuous-time random walk
        and correlated random walk state-space models to filter Argos
        satellite location data. Template Model Builder ('TMB') is
        used for fast estimation. The Argos data can be: (older) least
        squares-based locations; (newer) Kalman filter-based locations
        with error ellipse information; or a mixture of both. Separate
        measurement models are used for these two data types. The models
        estimate two sets of location states corresponding to: 1) each
        observation, which are (usually) irregularly timed; and 2)
        user-specified time intervals (regular or irregular).
    -   `r pkg("gazepath")` Eye-tracking data must be
        transformed into fixations and saccades before it can be
        analyzed. This package provides a non-parametric speed-based
        approach to do this on a trial basis. The method is especially
        useful when there are large differences in data quality, as the
        thresholds are adjusted accordingly. The same pre-processing
        procedure can be applied to all participants, while accounting
        for individual differences in data quality.
    -   `r pkg("marcher")` A set of tools for
        likelihood-based estimation, model selection and testing of two-
        and three-range shift and migration models for animal movement
        data as described in Gurarie et al. (2017). Provided movement
        data (X, Y and Time), including irregularly sampled data,
        functions estimate the time, duration and location of one or two
        range shifts, as well as the ranging area and auto-correlation
        structure of the movment. Tests assess, for example, whether the
        shift was "significant", and whether a two-shift migration was
        a true return migration.
    -   `r pkg("mdftracks")` 'MTrackJ' is an 'ImageJ'
        plugin for motion tracking and analysis. This package reads and
        writes 'MTrackJ Data Files' ('.mdf'). It supports 2D data
        and read/writes cluster, point, and channel information. If
        desired, generates track identifiers that are unique over the
        clusters. See the project page for more information and
        examples.
    -   `r pkg("mkde")` Provides functions to compute and
        visualize movement-based kernel density estimates (MKDEs) for
        animal utilization distributions in 2 or 3 spatial dimensions.
    -   `r pkg("momentuHMM")` Extended tools for analyzing
        telemetry data using generalized hidden Markov models. Features
        of momentuHMM (pronounced ``momentum'') include data
        pre-processing and visualization, fitting HMMs to location and
        auxiliary biotelemetry or environmental data, biased and
        correlated random walk movement models, multiple imputation for
        incorporating location measurement error and missing data,
        user-specified design matrices and constraints for covariate
        modelling of parameters, decoding of the state process,
        visualization of fitted models, model checking and selection,
        and simulation. See McClintock and Michelot (2018).
    -   `r pkg("mousetrack")` Extract from two-dimensional
        x-y coordinates of an arm-reaching trajectory, several dependent
        measures such as area under the curve, latency to start the
        movement, x-flips, etc.; which characterize the action-dynamics
        of the response. Mainly developed to analyze data coming from
        mouse-tracking experiments.
    -   `r pkg("mousetrap")` Mouse-tracking, the analysis of
        mouse movements in computerized experiments, is a method that is
        becoming increasingly popular in the cognitive sciences. The
        mousetrap package offers functions for importing, preprocessing,
        analyzing, aggregating, and visualizing mouse-tracking data.
    -   `r pkg("move")` Contains functions to access
        movement data stored in 'movebank.org' as well as tools to
        visualize and statistically analyze animal movement data, among
        others functions to calculate dynamic Brownian Bridge Movement
        Models. Move helps addressing movement ecology questions.
    -   `r pkg("movecost")` Provides the facility to
        calculate non-isotropic accumulated cost surface and least-cost
        paths using a number of human-movement-related cost functions
        that can be selected by the user. It just requires a Digital
        Terrain Model, a start location and (optionally) destination
        locations.
    -   `r pkg("moveHMM")` Provides tools for animal
        movement modelling using hidden Markov models. These include
        processing of tracking data, fitting hidden Markov models to
        movement data, visualization of data and fitted model, decoding
        of the state process.
    -   `r pkg("moveVis")` Tools to visualize movement data
        (e.g. from GPS tracking) and temporal changes of environmental
        data (e.g. from remote sensing) by creating video animations.
    -   `r pkg("moveWindSpeed")` Estimating wind speed from
        trajectories of individually tracked birds using a maximum
        likelihood approach.
    -   `r pkg("oce")` Supports the analysis of
        Oceanographic data, including 'ADCP' measurements,
        measurements made with 'argo' floats, 'CTD' measurements,
        sectional data, sea-level time series, coastline and topographic
        data, etc. Provides specialized functions for calculating
        seawater properties such as potential temperature in either the
        'UNESCO' or 'TEOS-10' equation of state. Produces graphical
        displays that conform to the conventions of the Oceanographic
        literature. This package is discussed extensively in Dan
        Kelley's book Oceanographic Analysis with R, published in 2018
        by 'Springer-Verlag' with ISBN 978-1-4939-8842-6.
    -   `r pkg("opentraj")` opentraj uses the Hybrid Single
        Particle Lagrangian Integrated Trajectory Model (HYSPLIT) for
        computing simple air parcel trajectories. The functions in this
        package allow users to run HYSPLIT for trajectory calculations,
        as well as get its results, directly from R without using any
        GUI interface.
    -   `r pkg("rerddapXtracto")` Contains three functions
        that access environmental data from any 'ERDDAP' data web
        service. The rxtracto() function extracts data along a
        trajectory for a given "radius" around the point. The
        rxtracto\_3D() function extracts data in a box. The
        rxtractogon() function extracts data in a polygon. All of those
        three function use the 'rerddap' package to extract the data,
        and should work with any 'ERDDAP' server. There are also two
        functions, plotBBox() and plotTrack() that use the 'plotdap'
        package to simplify the creation of maps of the data.
    -   `r pkg("riverdist")` Reads river network shape files
        and computes network distances. Also included are a variety of
        computation and graphical tools designed for fisheries telemetry
        research, such as minimum home range, kernel density estimation,
        and clustering analysis using empirical k-functions with a
        bootstrap envelope. Tools are also provided for editing the
        river networks, meaning there is no reliance on external
        software.
    -   `r pkg("SimilarityMeasures")` Functions to run and
        assist four different similarity measures. The similarity
        measures included are: longest common subsequence (LCSS),
        Frechet distance, edit distance and dynamic time warping (DTW).
        Each of these similarity measures can be calculated from two
        n-dimensional trajectories, both in matrix form.
    -   `r pkg("SiMRiv")` Provides functions to generate and
        analyze spatially-explicit individual-based multistate movements
        in rivers, heterogeneous and homogeneous spaces. This is done by
        incorporating landscape bias on local behaviour, based on
        resistance rasters. Although originally conceived and designed
        to simulate trajectories of species constrained to linear
        habitats/dendritic ecological networks (e.g. river networks),
        the simulation algorithm is built to be highly flexible and can
        be applied to any (aquatic, semi-aquatic or terrestrial)
        organism, independently on the landscape in which it moves.
        Thus, the user will be able to use the package to simulate
        movements either in homogeneous landscapes, heterogeneous
        landscapes (e.g. semi-aquatic animal moving mainly along rivers
        but also using the matrix), or even in highly contrasted
        landscapes (e.g. fish in a river network). The algorithm and its
        input parameters are the same for all cases, so that results are
        comparable. Simulated trajectories can then be used as
        mechanistic null models (Potts & Lewis 2014) to test a variety
        of 'Movement Ecology' hypotheses (Nathan et al. 2008),
        including landscape effects (e.g. resources, infrastructures) on
        animal movement and species site fidelity, or for predictive
        purposes (e.g. road mortality risk, dispersal/connectivity). The
        package should be relevant to explore a broad spectrum of
        ecological phenomena, such as those at the interface of animal
        behaviour, management, landscape and movement ecology, disease
        and invasive species spread, and population dynamics.
    -   `r pkg("smam")` Animal movement models including
        moving-resting process with embedded Brownian motion according
        to Yan et al. (2014), Pozdnyakov et al. (2017), Brownian motion
        with measurement error according to Pozdnyakov et al. (2014),
        and moving-resting-handling process with embedded Brownian
        motion, Pozdnyakov et al. (2018).
    -   `r pkg("SpaTimeClus")` Mixture model is used to
        achieve the clustering goal. Each component is itself a mixture
        model of polynomial autoregressive regressions whose the
        logistic weights consider the spatial and temporal information.
    -   `r pkg("stampr")` Perform spatial temporal analysis
        of moving polygons; a longstanding analysis problem in
        Geographic Information Systems. Facilitates directional
        analysis, shape analysis, and some other simple functionality
        for examining spatial-temporal patterns of moving polygons.
    -   `r pkg("surveillance")` Statistical methods for the
        modeling and monitoring of time series of counts, proportions
        and categorical data, as well as for the modeling of
        continuous-time point processes of epidemic phenomena. The
        monitoring methods focus on aberration detection in count data
        time series from public health surveillance of communicable
        diseases, but applications could just as well originate from
        environmetrics, reliability engineering, econometrics, or social
        sciences. The package implements many typical outbreak detection
        procedures such as the (improved) Farrington algorithm, or the
        negative binomial GLR-CUSUM method of Höhle and Paul (2008).
    -   `r pkg("trackdem")` Obtain population density and
        body size structure, using video material or image sequences as
        input. Functions assist in the creation of image sequences from
        videos, background detection and subtraction, particle
        identification and tracking. An artificial neural network can be
        trained for noise filtering. The goal is to supply accurate
        estimates of population size, structure and/or individual
        behavior, for use in evolutionary and ecological studies.
    -   `r pkg("trackdf")` Data frame class for storing
        collective movement data (e.g. fish schools, ungulate herds,
        baboon troops) collected from GPS trackers or computer vision
        tracking software.
    -   `r pkg("trackeR")` Provides infrastructure for
        handling running, cycling and swimming data from GPS-enabled
        tracking devices within R. The package provides methods to
        extract, clean and organise workout and competition data into
        session-based and unit-aware data objects of class
        'trackeRdata' (S3 class). The information can then be
        visualised, summarised, and analysed through flexible and
        extensible methods. Frick and Kosmidis (2017) \>, which is
        updated and maintained as one of the vignettes, provides
        detailed descriptions of the package and its methods, and
        real-data demonstrations of the package functionality.
    -   `r pkg("trackeRapp")` Provides an integrated user
        interface and workflow for the analysis of running, cycling and
        swimming data from GPS-enabled tracking devices.
    -   `r pkg("TrackReconstruction")` Reconstructs animal
        tracks from magnetometer, accelerometer, depth and optional
        speed data. Designed primarily using data from Wildlife
        Computers Daily Diary tags deployed on northern fur seals.
    -   `r pkg("trajectories")` Classes and methods for
        trajectory data, with support for nesting individual Track
        objects in track sets (Tracks) and track sets for different
        entities in collections of Tracks. Methods include selection,
        generalization, aggregation, intersection, simulation, and
        plotting.
    -   `r pkg("trajr")` A toolbox to assist with
        statistical analysis of 2-dimensional animal trajectories. It
        provides simple access to algorithms for calculating and
        assessing a variety of characteristics such as speed and
        acceleration, as well as multiple measures of straightness or
        tortuosity. McLean & Skowron Volponi (2018).
    -   `r pkg("trip")` Functions for accessing and
        manipulating spatial data for animal tracking, with
        straightforward coercion from and to other formats. Filter for
        speed and create time spent maps from animal track data. There
        are coercion methods to convert between 'trip' and 'ltraj'
        from 'adehabitatLT', and between 'trip' and 'psp' and
        'ppp' from 'spatstat'. Trip objects can be created from raw
        or grouped data frames, and from types in the 'sp', 'sf',
        'amt', 'trackeR', 'mousetrap', and other packages.
    -   `r pkg("tripEstimation")` Data handling and
        estimation functions for animal movement estimation from
        archival or satellite tags. Helper functions are included for
        making image summaries binned by time interval from Markov Chain
        Monte Carlo simulations.
    -   `r pkg("wildlifeDI")` Dynamic interaction refers to
        spatial-temporal associations in the movements of two (or more)
        animals. This package provides tools for calculating a suite of
        indices used for quantifying dynamic interaction with wildlife
        telemetry data. For more information on each of the methods
        employed see the references within. The package (as of version
        0.3) also has new tools for automating contact analysis in large
        tracking datasets. The package draws heavily on the classes and
        methods developed in the 'adehabitat' packages.

### Visualization

-   `r pkg("rasterVis")` includes a variety of methods that
    take advantage of the `z` slot of a `RasterStack` or `RasterBrick`
    object. Its [webpage](http://rastervis.r-forge.r-project.org/)
    includes several examples, from the hovmoller plot and horizon
    graph, to the density and histogram plots.
-   package `r pkg("googleVis")` provides an interface to
    show R data (tables) in the [Google Chart
    Tools](https://developers.google.com/chart/) .
    `r pkg("spacetime")` has a
    [vignette](https://cran.r-project.org/web/packages/spacetime/vignettes/stgvis.html)
    demonstrating its use for spatio-temporal data.
-   Package `r pkg("splancs")` provides animation and 3D
    interactive plots (using `r pkg("rgl")`) for displaying
    spatio-temporal point patterns.
-   `r pkg("mvtsplot")` provides multivariate time series
    plots, with examples on spatio-temporal data, published by [Peng
    (2008, Journal of Statistical
    Software)](http://www.jstatsoft.org/v25/c01/) .

### Data sets

-   Table data for fitting linear panel models are found in
    `r pkg("plm")`.
-   Package `r pkg("cshapes")` contains a data base with
    country boundaries, varying over time.
-   `r pkg("gstat")` contains the classic Irish wind data.
-   `r pkg("spacetime")` contains rural PM10 air quality
    measurements over Germany.
-   Some parts of the Cressie and Wikle (2011) book "Statistics for
    spatio-temporal data" can be reproduced by `demo(CressieWikle)` in
    `r pkg("spacetime")`.

### Retrieving data

Packages for retrieving data are:

-   Package `r pkg("openair")` has tools to analyze,
    interpret and understand air pollution data, but also tools to
    download UK air quality data.
-   `r pkg("ncdf4")` and `r pkg("RNetCDF")`
    allow reading and writing
    [netcdf](http://www.unidata.ucar.edu/software/netcdf/) files.
-   `r pkg("M3")` contains functions to read in and
    manipulate air quality model output from Models3-formatted files.
    This format is used by the Community Multiscale Air Quaility (CMAQ)
    model.
-   `r pkg("rmatio")` is a package for reading and writing
    Matlab MAT files from R.



### Links
-   [Pebesma (2012). spacetime: Spatio-Temporal Data in R. Journal of Statistical Software, 51(7).](https://www.jstatsoft.org/v51/i07/)
-   [Peng (2008). A Method for Visualizing Multivariate Time Series Data. Journal of Statistical Software, Code Snippets, 25(1).](https://www.jstatsoft.org/v25/c01/)
-   [Gapminder world: A tool for visualization of economic and health indicators per country (including access to the underlying tables in time-wide form).](http://www.gapminder.org/world/)
-   [Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019), Spatio-Temporal Statistics with R, Boca Raton, FL: Chapman & Hall/CRC (free to download).](https://spacetimewithr.org/)
