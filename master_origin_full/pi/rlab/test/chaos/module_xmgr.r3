//
// creating a xmgrace file for this
//
// setup page
    xmgfilename("gspc.gr",0);
    xmgpaper("letter", "portrait", 0);
//setup graph
    xmgtitle    ("Basic Time Series Analysis - S&P500 Market Index (gspc)", 1.3, 0);
    xmgsubtitle ("for period from 4-Dec-87 to 24-Mar-06", 1.3, 0);
    xmgxylabels (["", "Avg. Mutual Info."], 1.25, 0);
    xmgrange    ([ 0, NMAX, 0, 0.4 ], 0);
    xmgxticks   ([ 1000,9], ,0);
    xmgyticks   ([0.1,1],  ,0);
    xmggraphsize([0,1,0.75,1.0], 0);

    xmgxylabels (["Delay / Days", "Autocorrelation"], 1.25, 1);
    xmgrange    ([ 0, NMAX, -0.1, 0.8 ], 1);
    xmgxticks   ([ 1000,9], ,1);
    xmgyticks   ([0.1,1],  ,1);
    xmggraphsize( [0,1,0.4,0.65], 1);

    xmgxylabels (["Embedding Dimension", "False Nearest Neighbours"], 1.25, 2);
    xmgrange    ([ 1, 20, 0.0, 1.0 ], 2);
    xmgxticks   ([ 1,0], ,2);
    xmgyticks   ([0.2,1], ,2);
    xmggraphsize([0,1,0.0,0.25], 2);
    xmglegend  (  [0.48, 0.37], 1.0, 2);

// define the datasets
    xmgdataset  (y.[1], "xy", "change in closing price", 0);
    xmgstyle    ("line", [1,1,1.5], 0);
    xmgdataset  (y.[2], "xy", "volatility", 1);
    xmgstyle    ("line", [1,7,1.5], 1);
    xmgdataset  (y.[3], "xy", "change in closing price", 2);
    xmgstyle    ("line", [1,1,1.5], 2);
    xmgdataset  (y.[4], "xy", "volatility", 3);
    xmgstyle    ("line", [1,7,1.5], 3);
    xmgdataset  (y.[5], "xy", "change in closing price", 4);
    xmgstyle    ("line", [1,1,1.5], 4);
    xmgdataset  (y.[6], "xy", "volatility", 5);
    xmgstyle    ("line", [1,7,1.5], 5);
    xmggraph    ([0,1], 0);
    xmggraph    ([2,3], 1);
    xmggraph    ([4,5], 2);
    xmgpage     ([0,1,2],0);
    xmgprint    (0);



