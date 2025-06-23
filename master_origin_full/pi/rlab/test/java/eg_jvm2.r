//
//
//

jvm.create();

j_method = <<>>;
j_method.method     = "java/io/PrintStream/println";
j_method.descriptor = "(Ljava/lang/String;)V";

jvm.call_method(j_method, "from RLaB");


