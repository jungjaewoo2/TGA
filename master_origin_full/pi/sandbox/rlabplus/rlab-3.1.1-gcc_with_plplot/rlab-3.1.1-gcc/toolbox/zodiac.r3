//-------------------------------------------------------------------//
//
//  Synopsis:   Find the Chinese zodiac animal of a certain year
//
//  Syntax:     zodiac (year)
//
//  Example:
//
//  > rfile zodiac
//  > zodiac(1998)
//  Tiger  
//  > 
//
//  Note: Chinese new year day is usually in late January or early
//        February.  If you were born in January, your are probably 
//        the animal of the previous year.
//
//  T. S. Yang (yang@isec.com) 1/28/98  Happy New Year
//-------------------------------------------------------------------//
zodiac = function ( year )
{
   animal = ["Mouse", "Ox",   "Tiger",  "Rabbit",  "Dragon", "Snake",...   
             "Horse", "Goat", "Monkey", "Chicken", "Dog",    "Hog"]; 

   n = mod(year-1900, 12) + 1;

   return animal[n];
};
