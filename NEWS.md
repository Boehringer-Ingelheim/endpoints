# endpoints

## Development version

- `plot.makeDataSim()` no longer depends on `GGally` and now draws its plot
  matrix directly. As a result, it returns `invisible(NULL)` rather than a
  ggplot object, so plots can no longer be customized afterward with ggplot2's
  `+` syntax.
- `plot.makeDataSim()` gains a `title` argument for overriding the default
  `"Arm <k>"` plot title.
