using Pkg

Pkg.activate(".")
using JLD2
using GMT

ff =  jldopen("st_locs.jld2", "r")
lats = ff["lats"]
lons = ff["lons"]
sts = ff["st_idx"]
split("y 1 2 3 4 5 6 7 8 9", " ")
sts = "s0" .* split("2 3 1 4 10 5 7 6 9 8", " ")
sites = hcat(lons, lats)
close(ff)

main_region = (28.75, 29.75, -2.25, -1.25)
region = (28.75, 29.75, -2.25, -1.25)
G = GMT.gmtread("@earth_relief_01s", region=region)  # Example: NorCal
inset_region = (10, 50, -25, 15)
G_inset = GMT.gmtread("@earth_relief_03m", region= inset_region)  # Example: NorCal

# plot!(region= (28.25, 30.25, -3., -1.), pen = "1p,black", show = true)
rect_coords = [
            main_region[1] main_region[3];
            main_region[2] main_region[3];
            main_region[2] main_region[4];
            main_region[1] main_region[4];
            main_region[1] main_region[3]
]

inset_region = (10, 50, -25, 15)
cb = makecpt(cmap = :dem1, range = (500, 4500))
shade = grdgradient(G, azimuth=135, normalization=:t)

# inset

grdimage(G_inset, proj=:merc, cmap = cb, show=false, frame=:none)
coast!(region= inset_region, proj=:merc, shorelines="1/0.5p,black", ocean =:skyblue)
plot!(rect_coords, pen="3p,red", show = true) #, savefig = "inset_kivu_new.png")


gmtbegin()
# Overlay elevation with shading
GMT.grdimage(G, cmap= cb, show= false, shade= shade, frame = :auto)

GMT.colorbar!(position = (ouside =true, anchor =:MR), show = false, cmap = cb, frame=(annot=1000, label="Elevation [m]"))

coast!(borders="1/2p,black,--", show=false)  # Border thickness: 0.5pt
coast!(region=region, shorelines="1/0.5p,black",
      water=:skyblue, show=false)

text!(string.(sts), x=lons, y=lats, font=(10, "Helvetica", :black), offset=(-0., 0.25), legend = "MT sites", show = false)

feat_names= ["Rwanda", "Democratic", "Republic", "of Congo", "Uganda", "Burundi"]
feat_lats = 0.05 .+ [-2., -1.7, -1.75, -1.8, -1.08, -2.85]
feat_lons = [29.55, 28.85, 28.85, 28.85, 30, 29.75]
text!(feat_names, x=feat_lons, y=feat_lats, font=(14, "Helvetica-Bold", :black), offset=(-0., 0.25), show = false)

# volcanoes and lakes
feat_names= ["Mt Nyiragongo", "Mt Nyamuragira", "Lake Kivu"]
feat_lats = [-1.5, -1.4, -1.85]
feat_lons = [29.25, 29.2, 29.13]
text!(feat_names, x=feat_lons, y=feat_lats, font=(10, "Helvetica-Bold", :black), offset=(-0., 0.25), show = false)


inset(inset_box=(anchor=:TL, width=3, offset=(0.0, 0.0)), box=(fill=:white, pen=1))
    viz("inset_kivu_new.png", frame = :none, show = false)
inset(:end)

scatter!(lons, lats, marker = :i, markersize=0.3, mc=:black, markeredgecolor=:black, show=false) #, savefig = "kivu_sts_new.pdf")
gmtend(:show)