<?php
  include "corefunctions.php";
?>

<?php pageHeader("Documentation"); ?>

<body>

<div id="pagewrapper">

	<?php pageMenu("documentation"); ?>

	<div id="contentwrapper">
		<div id="leftside">
			<h3>PLplot Documentation</h3>
			<ul class="arrowlist">
	      <li> <a href="docbook-manual/">Browse the latest on-line documentation</a> </li>
	      <li> <a href="docbook-manual/.pdf">PDF form of documentation (<?php $size = number_format((filesize("docbook-manual/.pdf")/1024/1024),2); echo $size;?> MB)</a> </li>
	      <li> <a href="docbook-manual/.tar.gz">Info documentation tarball (<?php $size = number_format((filesize("docbook-manual/.tar.gz")/1024/1024),2); echo $size;?> MB)</a> </li>
	      <li> <a href="docbook-manual/.tar.gz">Man pages tarball (API chapter in nroff format) (<?php $size = number_format((filesize("docbook-manual/.tar.gz")/1024/1024),2); echo $size;?> MB)</a> </li>
	      <li> <a href="docbook-manual/.tar.gz">HTML results tarball (<?php $size = number_format((filesize("docbook-manual/.tar.gz")/1024/1024),2); echo $size;?> MB)</a> </li>
	      <li> <a href="doxygen/html">Browse doxygen-generated documentation</a> </li>
			</ul>
<p>
All but the doxygen-generated results above have been generated from
our DocBook/XML source files in doc/docbook/src using a variety of
DocBook backend tools.  The doxygen results above have been generated
by doxygen directly from a special form of comments in our source
code.  The DocBook and doxygen documentation builds happen
automatically (only on Linux since the required tools are only
available for that platform) with our CMake-based build system if you
specify the cmake options "-DBUILD_DOC=ON -DBUILD_DOX_DOC=ON".  For
more details about building the DocBook form of our documentation and
testing it, please look at the file doc/docbook/README.developers in
the source tree.
</p>
<p>
For those wishing to make some contribution to PLplot, helping out
with either/both the DocBook or doxygen documentation is a good place
to start.  For the DocBook case, the DocBook/XML syntax is quite
straightforward to understand if you simply follow the form of what is
already done in the files in the doc/docbook/src subdirectory of the
source tree. However, if you want to dig a little deeper into DocBook,
then this on-line book, <a
href="http://www.docbook.org/tdg/en/html/docbook.html">"DocBook: The
Definitive Guide"</a>, is an excellent reference.  For the doxygen
case, the documentation is controlled by a special form of comments in
our source files.  See src/pllegend.c for a good example of how to
document arguments of PLplot functions.  Note there are many other
files in src without this argument documentation at the present time
so there is plenty of scope to help out here.
</p>
		</div>

		<?php pageSidebar(1); ?>

		<div id="spacer"></div>
	</div>

	<?php pageFooter(); ?>
</div>

</body>
</html>
