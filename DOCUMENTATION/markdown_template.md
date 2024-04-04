# Template for Your Guide



# Chapter 1

Normal text, **bold text**, *italic text*

## Section 1

Bullet points: 
- one
- two
- three

Enumeration:
1. point
2. point
3. point
 

### How to display code

as multiline block (with copy to clipboard button)
```
sudo apt update
sudo apt install build-essential
sudo apt-get install manpages-dev
```
For multiline blocks you can also use syntax highlighting:
```python
def print_float(a : float):
    #indenting works just fine in the fenced code block
    s = f"float is {a:.2f}"
    print s
```

here is an example for bash syntax highlighting:

```bash
gmsh SphereWake . geo
```

Or as inline code block `module load gcc`

Or using an inline code block like that: 

    `module load gcc`

### how to add an image

place the file in the folder: `docs/source/_static`, to add an image use this code:

<img src="_static/BaseFlow1.png" alt="BaseFlow" width="400">

## how to do citations

The interested reader is referenced to Kaiser et al. [[1]](#1)

### References
<a id="1">[1]</a> 
Thomas Ludwig Kaiser, Thierry Poinsot, and Kilian Oberleithner. “Stability and Sensitivity Analysis of
Hydrodynamic Instabilities in Industrial Swirled Injection Systems”. In: Journal of Engineering for Gas
Turbines and Power 140.5 (Jan. 2018), p. 051506. issn: 0742-4795. doi: 10.1115/1.4038283. url: http:
//gasturbinespower.asmedigitalcollection.asme.org/article.aspx?doi=10.1115/1.4038283.

### Subcection 1

You can add links to [websites](https://stackoverflow.com), or to sections in this file [Link to Chapter 1](#chapter-1), or to other files [Link to README](README.md)

# Chapter 2

See [gitlab docs](https://docs.gitlab.com/ee/user/markdown.html) for more information on how to write markdown files

