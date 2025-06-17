//
//
//

jvm.create();

"\nTest 1:\n"
j_method = <<>>;
j_method.method     = "Prog.sum_me";
j_method.descriptor = "([D)D";

x = rand(1,6);
jvm.call_method(j_method, x)
sum(x)

"\nTest 2:\n"
j_method2 = <<>>;
j_method2.method     = "Prog.min_max";
j_method2.descriptor = "([D)[D";

x = rand(1,6);
jvm.call_method(j_method2, x)
[min(x),max(x)]

"\nTest 3:\n"
j_method3 = <<>>;
j_method3.method     = "Prog.min_len";
j_method3.descriptor = "([Ljava/lang/String;)Ljava/lang/String;";
x=["what is up?", "Not much.", "How about you?", "Help, I don't know what I am doing!"];
"the shortest string is " + jvm.call_method(j_method3, x)

"\nTest 4:\n"
j_method4 = <<>>;
j_method4.method     = "Prog.min_max_len";
j_method4.descriptor = "([Ljava/lang/String;)[Ljava/lang/String;";
x=["what is up?", "Not much.", "How about you?", "Help, I don't know what I am doing!"];
y = jvm.call_method(j_method4, x);
"the shortest string is " + y[1]
"the longest  string is " + y[2]

"\nTest 5: (This one should Fail - We are sending strings to function that expects doubles!)\n"
j_method5 = <<>>;
j_method5.method     = "Prog.min_max_len";
j_method5.descriptor = "([Ljava/lang/String;)[D";
x=["what is up?", "Not much.", "How about you?", "Help, I don't know what I am doing!"];
y = jvm.call_method(j_method5, x);
"the shortest string is " + y[1]
"the longest  string is " + y[2]
