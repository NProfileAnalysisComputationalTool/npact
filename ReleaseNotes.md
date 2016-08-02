# Version 3.0.0rc1 (zoomwindow) 2016-08-02
* Improved zoom window
 * improved gene display with codon highlighting
 * Gene navigation and editing
 * complementary strand bug fixes
 * saves work in progress via the url and historic track files
* Improved trac manipulation (drag and drop)
* UI polish and bug fixes

# Version 3.0pre
* Change backend to use Flask+GEvent
* Mycoplasma checkbox affects acgt_gamma and protein translation
* Use Bootstrap on frontend for styling


# Version 2.3.1 - 2015-07-30
* Fix direct strand protein translation

# Version 2.3.0 - 2015-07-30
* Blast a gene from the gene info window
* Scroll left and right with arrow keys
* Improve hit track rendering, show significance
* color blind friendly persists as a cookie
* Add modified ORFs track
* Draw ORF type tracks with varying styles to indicate type
* Goto base
* Search gene name
* Keep visible graph row at the top when the row height changes (adding/removing tracks)

# Version 2.2.1 - 2015-03-09
* polyfill for missing getChild in python2.6
* better logging of acgt_gamma exit
* Quiet some warnings in publish.sh

# Version 2.2.0 - 2015-03-04
* Ability to request PDF in new interface
* Print rendering of the web graphs
* config page is gone and that's all done in a panel on the results screen
* acgt_gamma and Allplots improvements to Hits drawing for significance

# Version 2.1.0

* graphs resize as the window resizes #843
* Render only visible graphs as you scroll #846
* toggle tracks visibility from the control panel #849
* clean the file (gbk,fasta) before doing anything #844
* render background shading in graphs #848
* can select which nucleotides to analyze #841
* bug fixes
   * correct scaling of nprofile lines
