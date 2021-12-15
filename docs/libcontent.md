# Provided Filter library

This page shows the content of the `pyphot` library with the respective
properties of the passband filters. The code to generate the table is also
provided in `hdf5test.cpp`.

## Table description

* `name`:  the identification name of the filter in the library.
* `detector type`: energy or photon counter.
* `wavelength units`:  filter defined with these units and all wavelength
  properties: `central wavelength`, `pivot wavelength`, and `effective wavelength`.
* `<X> mag`: magnitude in Vega, AB or ST system (w.r.t. the detector type)
* `<X> flux`: flux in :math:`erg/s/cm^2/\unicode{x212B}` in the `X` system
* `<X> Jy`: flux in :math:`Jy` (Jansky) in the `X` system


filters_properties.csv


@htmlonly
<div id="filterTable"></div>
<script src="https://d3js.org/d3.v3.min.js"></script>
<script type="text/javascript"charset="utf-8">
    d3.text("filters_properties.csv", function(data) {
        var parsedCSV = d3.csv.parseRows(data);

        var container = d3.select("contents")
            .append("table")

            .selectAll("tr")
                .data(parsedCSV).enter()
                .append("tr")

            .selectAll("td")
                .data(function(d) { return d; }).enter()
                .append("td")
                .text(function(d) { return d; });
    });
</script>
@endhtmlonly
