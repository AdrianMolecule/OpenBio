There are no Open Source Applications similar to Benchling, SnapGene or Geneious that work and I'm aware of.
This is an Open Source app to manipulate DNA and AA sequences.
It's based on a fairly rich library called BioPython.
It was started because Benchling is proprietary and does not allow extending the functionality.
I needed to qualitatively verify LAMP primers and there was no tool to do that.
At this point it acts as a simulator and can be use to simulate PCR, partial LAMP and I'll add Gibson assembly functionality soon. 
I added a lot of features like rotating the text for 3 to5 strands and show hydrogen bonds so it shoul be easy to use both in design and education.

There are some known limitations at this point. It supports just gb format but I also have unrelease code that do embl. Overlapping features don't visually overlap well, you can aneal one primer at a time and a strand with 2 loops is not very well depicted on one site.
Please contact me if you want me to improve on those.
The code was developed in VSS so the project should be pretty easy to install but to make it simpler I created a Windows Executable called main.exe you can download.
Once you click on the executable it takes several seconds to start but it should be fine afterwards.

The intention was to add features when I feel the need.

