// Translated by m2r alpha (version 0.52)
// Rlab target version 2.0.15
// # matlab filter(b,a,x)
// # rlab filter (  X , B, A )

runavg = function(x,w){
  f=filter(x,ones(1,w)/w,1).y; 
  // m2r// f=filter(ones(1,w)/w,1,x).y; 
  return f; 
  };
