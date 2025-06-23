//
// it all starts with hello world
//
public class Prog
{
  public static void gout0(String arg)
  {
    System.out.println("Hello World - Strings Again - gout0!");
    System.out.println("arg = " + arg);
  }
  public static void gout1(String[] args)
  {
    System.out.println("Hello World - Strings Again - gout1!");
    int i=0;
    for (String arg:args)
      System.out.println("args[" + (++i) + "] = " + arg);
  }
  public static void gout2(double[] vals)
  {
    int i=0;
    double result = 0;
    System.out.println("Hello World - Doubles!");
    for (double value:vals)
    {
      System.out.println("args[" + (++i) + "] = " + value);
      result += value;
    }
    System.out.println("sum is " + result);
  }
  public static double sum_me(double[] vals)
  {
    int i=0;
    double result = 0;
    System.out.println("Hello World : sum_me!");
    for (double value:vals)
    {
      System.out.println("sum_me: args[" + (++i) + "] = " + value);
      result += value;
    }
    System.out.println("sum_me: sum is " + result);
    return result;
  }
  public static double[] min_max(double[] vals)
  {
    System.out.println("Hello World : min_max!");
    int i=0;
    double[] rval;
    rval = new double[2];
    rval[0] = vals[0];
    rval[1] = vals[0];
    for (double value:vals)
    {
      System.out.println("min_max: args[" + (++i) + "] = " + value);
      if (rval[0]>value)
      {
        rval[0] = value;
        continue;
      }
      if (rval[1]<value)
      {
        rval[1] = value;
        continue;
      }
    }
    return rval;
  }
  public static String min_len(String[] args)
  {
    System.out.println("Hello World - min_len!");
    String rval;
    int i=0;
    rval = args[0];
    for (String value:args)
    {
      System.out.println("min_len: args[" + (++i) + "] = " + value);
      if (rval.length()>value.length())
      {
        rval = value;
        continue;
      }
    }
    System.out.println("min_len: the shortest string = " + rval);
    return rval;
  }
  public static String[] min_max_len(String[] args)
  {
    System.out.println("Hello World - min_max_len!");
    String[] rval;
    int i=0;
    rval = new String[2];
    rval[0] = args[0];
    rval[1] = args[0];
    for (String value:args)
    {
      System.out.println("min_max_len: args[" + (++i) + "] = " + value);
      if (rval[0].length()>value.length())
      {
        rval[0] = value;
        continue;
      }
      if (rval[1].length()<value.length())
      {
        rval[1] = value;
        continue;
      }
    }
    return rval;
  }
  public static void main(String[] args)
  {
    System.out.println("Hello World - Strings!");
    int i=0;
    for (String arg:args)
      System.out.println("args[" + (++i) + "] = " + arg);
  }
}

