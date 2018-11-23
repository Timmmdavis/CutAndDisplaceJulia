function triplot(vertices,faces,facecolor;opacity=0.5,cameraEye=[-1 -1 1])
    fcStr = [@sprintf("rgb(%5.2f,%5.2f,%5.2f)",255*facecolor[k,1],255*facecolor[k,2],255*facecolor[k,3]) for k = 1:size(facecolor,1) ]
    t = mesh3d(
        x = vertices[k,1], y = vertices[k,2], z = vertices[k,3],
        i = faces[k,1]-1,  j = faces[k,2]-1,  w = faces[k,3]-1, # Zero-based indexing in plotly
        facecolor=fcStr,opacity=opacity)
        layout = Layout(;title="Basic Triangle Plot",
                        scene=attr(;camera=attr(eye=attr(x=cameraEye[1],y=cameraEye[2],z=cameraEye[3]))))
    plot(t, layout)
end

##using Printf
##using PlotlyJS
##import PyPlot.triplot