//
//
//

jvm.create();

j_method0 = <<>>;
j_method0.method     = "Prog.gout0";
j_method0.descriptor = "(Ljava/lang/String;)V";

jvm.call_method(j_method0, ["1: from RLaB"]);

j_method1 = <<>>;
j_method1.method     = "Prog.main";
j_method1.descriptor = "([Ljava/lang/String;)V";

jvm.call_method(j_method1, ["1: from RLaB"]);
jvm.call_method(j_method1, ["2: from RLaB", "2: more from RLaB", "2: even more from RLaB"]);

j_method2 = <<>>;
j_method2.method     = "Prog.gout1";
j_method2.descriptor = "([Ljava/lang/String;)V";

jvm.call_method(j_method2, ["1: from RLaB"]);
jvm.call_method(j_method2, ["2: from RLaB", "2: more from RLaB", "2: even more from RLaB"]);

j_method3 = <<>>;
j_method3.method     = "Prog.gout2";
j_method3.descriptor = "([D)V";

jvm.call_method(j_method3, rand(1,2));
jvm.call_method(j_method3, rand());
