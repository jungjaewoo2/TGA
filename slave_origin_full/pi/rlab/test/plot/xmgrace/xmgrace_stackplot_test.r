// test file for xmgrace.r toolkit

x=rand(1000,2);


// setup page
    xmgfilename("test.gr",0);
    xmgpaper("letter", "landscape", 0);
//setup graph
    xmgtitle    ("Fig.1: Plot of randomly generated points", 1.3, 0);
    xmgxylabels (["x", "y"], 0);
    xmgrange    ([ 0, 1, 0, 1 ], 0);
    xmgxticks   ([0.1,1], 0);
    xmgyticks   ([0.1,1],  0);
    xmggraphsize( [-0.1,0.5,-.1, 0.4], 0);
    xmggraphsize( [ 0.25,0.49,0.1,0.39], 1);
    xmggraphsize( [0.5,0.7,-0.1, 0.4], 2);
    xmggraphsize( [0,1,0.5, 0.8], 3);
    xmggraphsize( [0.6,0.9,0.8, 1.1], 4);
// define the datasets
    xmgdataset  ( x, "xy", "random points", 0);
    xmgstyle    ("symbol", [1,1,0.1], 0);
    xmgdataset  ( x, "xy", "random points", 1);
    xmgstyle    ("symbol", [1,2,0.1], 1);
    xmgdataset  ( x, "xy", "random points", 2);
    xmgstyle    ("symbol", [1,3,0.1], 2);
    xmgdataset  ( x, "xy", "random lines", 3);
    xmgstyle    ("line", [1,4,0.1], 3);
    xmgdataset  ( x, "xy", "random points", 4);
    xmgstyle    ("symbol", [1,5,0.1], 4);

    xmggraph    ([0], 0);
    xmggraph    ([1], 1);
    xmggraph    ([2], 2);
    xmggraph    ([3], 3);
    xmggraph    ([4], 4);
    xmgpage     ( [0:4],0);
    xmgprint    (0);
