* [Load the ggplot2 package](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#load-the-ggplot2-package)
* [Initialize a plot (`ggplot()`)](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#initialize-a-plot-ggplot)
* [Aesthetic mappings (`aes()`)](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#aesthetic-mappings-aes)
* [Adding components](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#adding-components)
* [Defaults and overrides](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#defaults-and-overrides)
* [Facets (`facet_wrap()` and `facet_grid()`)](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#facets-facet_wrap-and-facet_grid)
* [Changing cosmetic details](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#changing-cosmetic-details)
* [More](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#more)
* [Extending ggplot2](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#extending-ggplot2)
* [Decode](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#decode)
* [References](https://github.com/nmnh-r-users/meetups/tree/master/code/ggplot#references)

ggplot2 is a widely used package for data visualization, allowing you to create highly customizable plots with aesthetically pleasing defaults.  Documentation for ggplot2 is quite good and many resources exist out there that can provide you with example code for various modifications.

The syntax for ggplot2 is different than what we've encountered before.  Think of it, however, as functioning much in the same way as Adobe Illustrator - you build your plots layer by layer, adjusting either each layer or the entire workspace until the final image looks the way you want it to.

Check out [Cookbook for R](http://www.cookbook-r.com/Graphs/) for worked examples and explanations of how to make minor modifications to plots.

# Load the ggplot2 package

Remember how?

# Initialize a plot (`ggplot()`)

First, initialize a new plot with the function `ggplot()`.  Often, this function is used to define shared defaults for all the other layers - for example, what dataset to reference or what should be used as your x and y coordinates.  We will get to this in the next step.

```R
ggplot()

# let's use the iris dataset for our plot
ggplot(data = iris)
```

Just like numbers and data frames, this plot is an object that we can save to a variable to recall later.

``` R
g <- ggplot(data = iris)

g

class(g)
```

# Aesthetic mappings (`aes()`)

The function `aes()` defines the link between data and visuals, controlling data-dependent aspects of the plot.  Think things like:
* x and y axis scales (should x go from -2 to 10 or from 0 to 500?)
* positions of points (what coordinates should points be drawn at?)

In these cases, the plot must refer to your data in order to learn how items should be displayed.

``` R
ggplot(data = iris)

# we'll add a default location for the plot to find x and y values
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width))
```

# Adding components

It is impossible to fully outline the purpose of `aes()` without first getting the hang of how one builds up components of a plot layer by layer.

We add layers using the plus sign (`+`).  There's several kinds of objects you can add to a plot:
* `geom_xxx()`
* `stat_xxx()`
* layers that adjust the overall look (kind of like settings or preferences)

Full list of geoms: http://ggplot2.tidyverse.org/reference/#section-layer-geoms

Let's add a layer of points:

``` R
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point()
```
And also a linear regression line:

```R
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point() +
  geom_smooth(method = lm)
```

Or maybe you want to use the x and y coordinates to plot text:

```R
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_text(aes(label = Species))
```

Notice each layer can have its own `aes()` function to define aspects that are unique to just that layer.

### Order matters

ggplot2 will render the plot as it reads it.  That means each added layer gets placed on top of the previous one.

>### Exercise
>
> Switch the order of `geom_point()` and `geom_smooth()`.  Notice the difference?

### To `aes()` or not to `aes()`

Visual features that do not depend on the data should be listed outside of `aes()`.  For example, if points in our iris plot should be colored by Species, then we add the argument "color" inside `aes()`:

``` R
ggplot(data = iris) +
  geom_point(aes(x = Sepal.Length, y = Sepal.Width, color = Species))
```

However, if **all** points should be magenta, then "color" should be outside of `aes()`:

``` R
ggplot(data = iris) +
  geom_point(aes(x = Sepal.Length, y = Sepal.Width), color = "magenta")
```

# Defaults and overrides

Placement and order are key.  I mentioned earlier that arguments in the initial function `ggplot()` are used, if applicable, as default for all other layers in the plot.

> ### Exercise
> 
> Turn to the person next to you and chat about:
> 1) what is different between the code in the following four plots
> 2) what effect do the changes have on how the plot is rendered (use the little arrows in the plot window of RStudio to flip back and forth between the plots)
> 3) why is it doing these things?
> 
> ``` R
> plot1 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
>   geom_point() +
>   geom_smooth(method = lm)
> 
> ###
> 
> plot2 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
>   geom_point(aes(colour = Species)) +
>   geom_smooth(method = lm)
> 
> ###
> 
> plot3 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, colour = Species)) +
>   geom_point() +
>   geom_smooth(method = lm)
> 
> ###
> 
> plot4 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, colour = Species)) +
>   geom_point(colour = "peachpuff3") +
>   geom_smooth(method = lm)
> ```

# Facets (`facet_wrap()` and `facet_grid()`)

If you need more convincing that ggplot2 is a powerful graphics tool and that you should use it, here it is.  Facets (facet wrap and facet grid) allow for easy organizing of plots.

Take a look at the dataset `CO2`.

```R
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant))

# compare the two types
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant, shape = Type)) +
  geom_smooth(aes(lty = Type))
```

That's tough to interpret.  `facet_wrap()` explodes the original plot out "as a function of (~)" some variable.  Each sub-plot is then presented one after the other.

```R
# here, the facet wrap variable is Plant
# you can see that every color (representing each Plant) has it's own sub-plot
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant)) +
  facet_wrap(~ Plant)

# here, the facet wrap variable is Type
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant)) +
  facet_wrap(~ Type)
```

`facet_grid()` does the same kind of subsetting that `facet_wrap()` does, except it makes a table using a combination of variables (2 dimensions instead of 1).

```R
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant)) +
  facet_grid(Treatment ~ Type)
```

# Changing cosmetic details

Every aspect of the plot you see can be turned on or off and adjusted.  We'll cover just a few major modification you can make.

http://ggplot2.tidyverse.org/articles/ggplot2-specs.html

### Themes

"Theme" refers broadly to cosmetic details - colors, font size, spacing, tick marks, etc.

There are a number of built-in pre-packaged themes that are quite nice looking without much modification.  Examples of those can be found [here](http://ggplot2.tidyverse.org/reference/ggtheme.html).  You add them as layers as you would geoms and stats.

```R
plot3 +
  theme_minimal()
```

The function `theme()` allows you to really fine-tune things.  It can take a lot of arguments, which are listed [here](http://ggplot2.tidyverse.org/reference/theme.html).  Usually, the arguments take a special function called `element_xxx()` that describes what the specific shape or text looks like.

For example, the background of a plot is in reality just a drawn rectangle with a fill of "white" and no observable border.  By defining the background as a rectangle with different properties (using the function `element_rect()`), we change the look of the plot.

```R
plot3 +
  theme_bw() # use the black and white theme

plot3 +
  theme_bw() +
  theme(axis.text.y = element_text(size = 38, colour = "blue"),
        plot.background = element_rect(fill = "violet"))
```

You can also set an argument as `element_blank()` to remove it entirely.

```R
# remove tickmark text from both axis
plot3 +
  theme_bw() +
  theme(axis.text = element_blank(),
        plot.background = element_rect(fill = "violet"))
```

> ### Exercise
> 
> 1. You see that grid of white lines in the background of plot3?  Change that grid to any color you'd like.
> 2. Re-position the legend in plot3 from the right hand side to the bottom.

### Scales

`scale_xxx_yyy()` functions adjust anything having to do with data boundaries: how far apart tick marks are, the range colors to use in a heatmap, the colors of categories, the minimum and maximum values to show on an axis, etc.

```R
# only show conc from 0 to 300
ggplot(data = CO2, aes(x = conc, y = uptake)) +
  geom_point(aes(colour = Plant)) +
  facet_grid(Treatment ~ Type) +
  scale_x_continuous(limit = c(0, 300))
```

> ### Exercise
> 
> 1. Reverse the x axis on plot3 from earlier.
> 
> 2. Change the three colors used in plot3 to "olivedrab", "sienna", and "lawngreen".

# More

You can plot multiple datasets on the same plot by defining the "data = " argument for each layer that does not match the default.

```R
text.corners <- data.frame(
  x = c(4, 4, 8, 6),
  y = c(2, 4, 2, 3),
  text = c("bottom-left", "bottom-right", "top-left", "center")
)

ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point() +
  geom_text(data = text.corners, aes(x = x, y = y, label = text))
```

> ### Exercise
> 
> Add the text "top-right" to the top right of the plot above.

# Extending ggplot2

There exist many packages that build on ggplot2 to provide shortcuts for visualizing and modifying specific content. 
 Remember, though, packages are just collections of functions.  All of these packages USE ggplot2 functions and the structure of ggplot2 objects.  So if you know ggplot2 and how flexible it is, you can update and customize any canned function.  Never be just content with the defaults.

Packages of interest (this list can grow):
* [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) - for visualizing phylogenetic trees, including annotations, figures, etc.
* [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html), [patchwork](https://github.com/thomasp85/patchwork), and others - arrange ggplot objects with panel labels and flexible spacing
* [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html) - offset labels to improve readability

# Decode

* What geoms are in this plot?
* How many datasets are being used?
* What is the order in which you have to add elements?

<img src="https://github.com/nmnh-r-users/R-Basics/blob/master/plot.div.png" width = "80%">

[annotated version](https://github.com/nmnh-r-users/R-Basics/blob/master/plot.div.annotate.png)

# References

[Cookbook for R's ggplot2 resource](http://www.cookbook-r.com/Graphs/)

[Reference page for ggplot2 functions](http://ggplot2.tidyverse.org/reference/)

[Cheat sheet](https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf)

\*\*[Breaking down a plot from the Economist](http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html#putting_it_all_together)\*\* (this is a fun one)
