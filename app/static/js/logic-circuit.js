var red = "orangered";  // 0 or false
var green = "forestgreen";  // 1 or true

function init() {
    var $ = go.GraphObject.make;  // for conciseness in defining templates

    myDiagram =
    $(go.Diagram, "myDiagramDiv",  // create a new Diagram in the HTML DIV element "myDiagramDiv"
        {
        "draggingTool.isGridSnapEnabled": true,  // dragged nodes will snap to a grid of 10x10 cells
        "undoManager.isEnabled": true
        });

    // when the document is modified, add a "*" to the title and enable the "Save" button
    myDiagram.addDiagramListener("Modified", function(e) {
    var button = document.getElementById("saveModel");
    if (button) button.disabled = !myDiagram.isModified;
    var idx = document.title.indexOf("*");
    if (myDiagram.isModified) {
        if (idx < 0) document.title += "*";
    } else {
        if (idx >= 0) document.title = document.title.substr(0, idx);
    }
    });

    var palette = new go.Palette("palette");  // create a new Palette in the HTML DIV element "palette"

    // creates relinkable Links that will avoid crossing Nodes when possible and will jump over other Links in their paths
    myDiagram.linkTemplate =
    $(go.Link,
        {
        routing: go.Link.AvoidsNodes,
        curve: go.Link.JumpOver,
        corner: 3,
        relinkableFrom: true, relinkableTo: true,
        selectionAdorned: false, // Links are not adorned when selected so that their color remains visible.
        shadowOffset: new go.Point(0, 0), shadowBlur: 5, shadowColor: "blue",
        },
        new go.Binding("isShadowed", "isSelected").ofObject(),
        $(go.Shape,
        { name: "SHAPE", strokeWidth: 2, stroke: red }));

    // node template helpers
    var sharedToolTip =
    $("ToolTip",
        { "Border.figure": "RoundedRectangle" },
        $(go.TextBlock, { margin: 2 },
        new go.Binding("text", "", function(d) { return d.category; })));

    // define some common property settings
    function nodeStyle() {
    return [new go.Binding("location", "loc", go.Point.parse).makeTwoWay(go.Point.stringify),
    new go.Binding("isShadowed", "isSelected").ofObject(),
    {
        selectionAdorned: false,
        shadowOffset: new go.Point(0, 0),
        shadowBlur: 15,
        shadowColor: "blue",
        toolTip: sharedToolTip
    }];
    }

    function shapeStyle() {
    return {
        name: "NODESHAPE",
        fill: "lightgray",
        stroke: "darkslategray",
        desiredSize: new go.Size(40, 40),
        strokeWidth: 2
    };
    }

    function portStyle(input) {
    return {
        desiredSize: new go.Size(6, 6),
        fill: "black",
        fromSpot: go.Spot.Right,
        fromLinkable: !input,
        toSpot: go.Spot.Left,
        toLinkable: input,
        toMaxLinks: 1,
        cursor: "pointer"
    };
    }

    // define templates for each type of node
    var inputTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "Circle", shapeStyle(),
        { fill: green }),  // override the default fill (from shapeStyle()) to be red
        $(go.Shape, "Rectangle", portStyle(false),  // the only port
        { portId: "", alignment: new go.Spot(1, 0.5) }),
        { // if double-clicked, an input node will change its value, represented by the color.
        doubleClick: function(e, obj) {
            e.diagram.startTransaction("Toggle Input");
            var shp = obj.findObject("NODESHAPE");
            shp.fill = (shp.fill === green) ? red : green;
            updateStates();
            e.diagram.commitTransaction("Toggle Input");
        }
        }
    );

    var outputTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "Rectangle", shapeStyle(),
        { fill: green }),  // override the default fill (from shapeStyle()) to be green
        $(go.Shape, "Rectangle", portStyle(true),  // the only port
        { portId: "", alignment: new go.Spot(0, 0.5) })
    );

    var andTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "AndGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var orTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "OrGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0.16, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0.16, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var xorTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "XorGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0.26, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0.26, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var norTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "NorGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0.16, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0.16, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var xnorTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "XnorGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0.26, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0.26, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var nandTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "NandGate", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in1", alignment: new go.Spot(0, 0.3) }),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in2", alignment: new go.Spot(0, 0.7) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    var notTemplate =
    $(go.Node, "Spot", nodeStyle(),
        $(go.Shape, "Inverter", shapeStyle()),
        $(go.Shape, "Rectangle", portStyle(true),
        { portId: "in", alignment: new go.Spot(0, 0.5) }),
        $(go.Shape, "Rectangle", portStyle(false),
        { portId: "out", alignment: new go.Spot(1, 0.5) })
    );

    // add the templates created above to myDiagram and palette
    myDiagram.nodeTemplateMap.add("input", inputTemplate);
    myDiagram.nodeTemplateMap.add("output", outputTemplate);
    myDiagram.nodeTemplateMap.add("and", andTemplate);
    myDiagram.nodeTemplateMap.add("or", orTemplate);
    myDiagram.nodeTemplateMap.add("xor", xorTemplate);
    myDiagram.nodeTemplateMap.add("not", notTemplate);
    myDiagram.nodeTemplateMap.add("nand", nandTemplate);
    myDiagram.nodeTemplateMap.add("nor", norTemplate);
    myDiagram.nodeTemplateMap.add("xnor", xnorTemplate);

    // share the template map with the Palette
    palette.nodeTemplateMap = myDiagram.nodeTemplateMap;

    palette.model.nodeDataArray = [
    { category: "input" },
    { category: "output" },
    { category: "and" },
    { category: "or" },
    { category: "xor" },
    { category: "not" },
    { category: "nand" },
    { category: "nor" },
    { category: "xnor" }
    ];

    // load the initial diagram
    load();

    // continually update the diagram
    loop();
}
// update the diagram every 250 milliseconds
function loop() {
    setTimeout(function() { updateStates(); loop(); }, 250);
}
// update the value and appearance of each node according to its type and input values
function updateStates() {
    var oldskip = myDiagram.skipsUndoManager;
    myDiagram.skipsUndoManager = true;
    // do all "input" nodes first
    myDiagram.nodes.each(function(node) {
    if (node.category === "input") {
        doInput(node);
    }
    });
    // now we can do all other kinds of nodes
    myDiagram.nodes.each(function(node) {
    switch (node.category) {
        case "and": doAnd(node); break;
        case "or": doOr(node); break;
        case "xor": doXor(node); break;
        case "not": doNot(node); break;
        case "nand": doNand(node); break;
        case "nor": doNor(node); break;
        case "xnor": doXnor(node); break;
        case "output": doOutput(node); break;
        case "input": break;  // doInput already called, above
    }
    });
    myDiagram.skipsUndoManager = oldskip;
}
// helper predicate
function linkIsTrue(link) {  // assume the given Link has a Shape named "SHAPE"
    return link.findObject("SHAPE").stroke === green;
}
// helper function for propagating results
function setOutputLinks(node, color) {
    node.findLinksOutOf().each(function(link) { link.findObject("SHAPE").stroke = color; });
}
// update nodes by the specific function for its type
// determine the color of links coming out of this node based on those coming in and node type
function doInput(node) {
    // the output is just the node's Shape.fill
    setOutputLinks(node, node.findObject("NODESHAPE").fill);
}
function doAnd(node) {
    var color = node.findLinksInto().all(linkIsTrue) ? green : red;
    setOutputLinks(node, color);
}
function doNand(node) {
    var color = !node.findLinksInto().all(linkIsTrue) ? green : red;
    setOutputLinks(node, color);
}
function doNot(node) {
    var color = !node.findLinksInto().all(linkIsTrue) ? green : red;
    setOutputLinks(node, color);
}
function doOr(node) {
    var color = node.findLinksInto().any(linkIsTrue) ? green : red;
    setOutputLinks(node, color);
}
function doNor(node) {
    var color = !node.findLinksInto().any(linkIsTrue) ? green : red;
    setOutputLinks(node, color);
}
function doXor(node) {
    var truecount = 0;
    node.findLinksInto().each(function(link) { if (linkIsTrue(link)) truecount++; });
    var color = truecount % 2 !== 0 ? green : red;
    setOutputLinks(node, color);
}
function doXnor(node) {
    var truecount = 0;
    node.findLinksInto().each(function(link) { if (linkIsTrue(link)) truecount++; });
    var color = truecount % 2 === 0 ? green : red;
    setOutputLinks(node, color);
}
function doOutput(node) {
    // assume there is just one input link
    // we just need to update the node's Shape.fill
    node.linksConnected.each(function(link) { node.findObject("NODESHAPE").fill = link.findObject("SHAPE").stroke; });
}
// save a model to and load a model from JSON text, displayed below the Diagram
function save() {
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
    myDiagram.isModified = false;
}
function load() {
    myDiagram.model = go.Model.fromJson(document.getElementById("mySavedModel").value);
}

////////// custom functions

function reset_design() {
    var init_design = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [],
        "linkDataArray": []
    }
    myDiagram.model = go.Model.fromJson(init_design);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();

    var design = document.getElementById("design_id");
    design.options[0].selected = true;
}

function load_model() {
    var model = document.getElementById("model_id");
    var selected = model.options[model.selectedIndex].value;
    var upload_model = document.getElementById('upload_model');
    if (selected=='Upload custom model') {
        upload_model.style.display = 'block';
    }
    else {
        upload_model.style.display = 'none';
    }
}

function load_builtin() {
    var design = document.getElementById("design_id");
    var selected = design.options[design.selectedIndex].value;
    var upload_circuit = document.getElementById('upload_circuit');
    if (selected=='Upload custom design') {
        upload_circuit.style.display = 'block';
        var init_design = {
            "class": "GraphLinksModel",
            "linkFromPortIdProperty": "fromPort",
            "linkToPortIdProperty": "toPort",
            "nodeDataArray": [],
            "linkDataArray": []
        }
        myDiagram.model = go.Model.fromJson(init_design);
        document.getElementById("mySavedModel").value = myDiagram.model.toJson();
    }
    else {
        upload_circuit.style.display = 'none';
        if (selected=='3-inputs AND gate') { load_and3(); }
        else if (selected=='4-inputs AND gate') { load_and4(); }
        else if (selected=='8-inputs AND gate') { load_and8(); }
        else if (selected=='1-bit half adder') { load_halfadder(); }
        else if (selected=='1-bit full adder') { load_fulladder(); }
    }
}

function load_and3() { //incorrect
    var and3 = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [
            {"category":"input","key":"input1","loc":"-150 -330"},
            {"category":"input","key":"input2","loc":"-150 -270"},
            {"category":"input","key":"input3","loc":"-150 -210"},
            {"category":"and","key":"and1","loc":"-50 -300"},
            {"category":"and","key":"and2","loc":"60 -240"},
            {"category":"output","key":"output1","loc":"170 -240"}
        ],
        "linkDataArray": [
            {"from":"input1","to":"and1","fromPort":"","toPort":"in1"},
            {"from":"input2","to":"and1","fromPort":"","toPort":"in2"},
            {"from":"and1","to":"and2","fromPort":"","toPort":"in1"},
            {"from":"input3","to":"and2","fromPort":"","toPort":"in2"},
            {"from":"and2","to":"output1","fromPort":"out","toPort":""}
        ]
    }
    myDiagram.model = go.Model.fromJson(and3);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
}

function load_and4() { //incorrect
    var and4 = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [
            {"category":"input","key":"input1","loc":"-150 -150"},
            {"category":"input","key":"input12","loc":"-150 -210"},
            {"category":"output","key":"output1","loc":"170 -240"},
            {"category":"input","key":-1,"loc":"-150 -270"},
            {"category":"input","key":-10,"loc":"-150 -330"},
            {"category":"and","key":-3,"loc":"-50 -180"},
            {"category":"and","key":-11,"loc":"-50 -300"},
            {"category":"and","key":-14,"loc":"60 -240"}
        ],
        "linkDataArray": [
            {"from":"input1","to":-3,"fromPort":"","toPort":"in2"},
            {"from":-1,"to":-11,"fromPort":"","toPort":"in2"},
            {"from":-10,"to":-11,"fromPort":"","toPort":"in1"},
            {"from":"input12","to":-3,"fromPort":"","toPort":"in1"},
            {"from":-6,"to":-12,"fromPort":"","toPort":"in1"},
            {"from":-7,"to":-12,"fromPort":"","toPort":"in2"},
            {"from":-8,"to":-13,"fromPort":"","toPort":"in1"},
            {"from":-9,"to":-13,"fromPort":"","toPort":"in2"},
            {"from":-11,"to":-14,"fromPort":"out","toPort":"in1"},
            {"from":-3,"to":-14,"fromPort":"out","toPort":"in2"},
            {"from":-12,"to":-15,"fromPort":"out","toPort":"in1"},
            {"from":-13,"to":-15,"fromPort":"out","toPort":"in2"},
            {"from":-14,"to":-16,"fromPort":"out","toPort":"in1"},
            {"from":-15,"to":-16,"fromPort":"out","toPort":"in2"},
            {"from":-14,"to":"output1","fromPort":"out","toPort":""}
        ]
    }
    myDiagram.model = go.Model.fromJson(and4);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
}

function load_and8() {
    var and8 = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [
            {"category":"input","key":"input1","loc":"-150 -150"},
            {"category":"output","key":"output1","loc":"280 -120"},
            {"category":"and","key":-3,"loc":"-50 -180"},
            {"category":"input","key":"input12","loc":"-150 -210"},
            {"category":"input","key":-1,"loc":"-150 -270"},
            {"category":"input","key":-6,"loc":"-150 -90"},
            {"category":"input","key":-7,"loc":"-150 -30"},
            {"category":"input","key":-8,"loc":"-150 30"},
            {"category":"input","key":-9,"loc":"-150 90"},
            {"category":"input","key":-10,"loc":"-150 -330"},
            {"category":"and","key":-11,"loc":"-50 -300"},
            {"category":"and","key":-12,"loc":"-50 -60"},
            {"category":"and","key":-13,"loc":"-50 60"},
            {"category":"and","key":-14,"loc":"60 -240"},
            {"category":"and","key":-15,"loc":"60 0"},
            {"category":"and","key":-16,"loc":"170 -120"}
        ],
        "linkDataArray": [
            {"from":"input1","to":-3,"fromPort":"","toPort":"in2"},
            {"from":-1,"to":-11,"fromPort":"","toPort":"in2"},
            {"from":-10,"to":-11,"fromPort":"","toPort":"in1"},
            {"from":"input12","to":-3,"fromPort":"","toPort":"in1"},
            {"from":-6,"to":-12,"fromPort":"","toPort":"in1"},
            {"from":-7,"to":-12,"fromPort":"","toPort":"in2"},
            {"from":-8,"to":-13,"fromPort":"","toPort":"in1"},
            {"from":-9,"to":-13,"fromPort":"","toPort":"in2"},
            {"from":-11,"to":-14,"fromPort":"out","toPort":"in1"},
            {"from":-3,"to":-14,"fromPort":"out","toPort":"in2"},
            {"from":-12,"to":-15,"fromPort":"out","toPort":"in1"},
            {"from":-13,"to":-15,"fromPort":"out","toPort":"in2"},
            {"from":-14,"to":-16,"fromPort":"out","toPort":"in1"},
            {"from":-15,"to":-16,"fromPort":"out","toPort":"in2"},
            {"from":-16,"to":"output1","fromPort":"out","toPort":""}
        ]
    }
    myDiagram.model = go.Model.fromJson(and8);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
}

function load_halfadder() {
    var half_adder = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [
            {"category":"output","key":"output1","loc":"150 -140"},
            {"category":"input","key":-1,"loc":"-140 -170"},
            {"category":"input","key":-10,"loc":"-140 -230"},
            {"category":"and","key":-11,"loc":"-40 -200"},
            {"category":"xor","key":-5,"loc":"40 -140"},
            {"category":"output","key":-2,"loc":"150 -200"}
        ],
        "linkDataArray": [
            {"from":-1,"to":-11,"fromPort":"","toPort":"in2"},
            {"from":-10,"to":-11,"fromPort":"","toPort":"in1"},
            {"from":-10,"to":-5,"fromPort":"","toPort":"in1"},
            {"from":-1,"to":-5,"fromPort":"","toPort":"in2"},
            {"from":-5,"to":"output1","fromPort":"out","toPort":""},
            {"from":-11,"to":-2,"fromPort":"out","toPort":""}
        ]
    }
    myDiagram.model = go.Model.fromJson(half_adder);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
}

function load_fulladder() {
    var full_adder = {
        "class": "GraphLinksModel",
        "linkFromPortIdProperty": "fromPort",
        "linkToPortIdProperty": "toPort",
        "nodeDataArray": [
            {"category":"input","key":-1,"loc":"-180 -240"},
            {"category":"input","key":-10,"loc":"-180 -300"},
            {"category":"and","key":-11,"loc":"-80 -270"},
            {"category":"xor","key":-5,"loc":"-80 -210"},
            {"category":"input","key":-6,"loc":"-180 -120"},
            {"category":"and","key":-3,"loc":"20 -150"},
            {"category":"xor","key":-7,"loc":"120 -100"},
            {"category":"or","key":-4,"loc":"120 -190"},
            {"category":"output","key":-2,"loc":"230 -190"},
            {"category":"output","key":-12,"loc":"230 -100"}
        ],
        "linkDataArray": [
            {"from":-10,"to":-11,"fromPort":"","toPort":"in1"},
            {"from":-1,"to":-11,"fromPort":"","toPort":"in2"},
            {"from":-10,"to":-5,"fromPort":"","toPort":"in1"},
            {"from":-1,"to":-5,"fromPort":"","toPort":"in2"},
            {"from":-11,"to":-4,"fromPort":"out","toPort":"in1"},
            {"from":-3,"to":-4,"fromPort":"out","toPort":"in2"},
            {"from":-5,"to":-3,"fromPort":"out","toPort":"in1"},
            {"from":-5,"to":-7,"fromPort":"out","toPort":"in1"},
            {"from":-6,"to":-3,"fromPort":"","toPort":"in2"},
            {"from":-6,"to":-7,"fromPort":"","toPort":"in2"},
            {"from":-4,"to":-2,"fromPort":"out","toPort":""},
            {"from":-7,"to":-12,"fromPort":"out","toPort":""}
        ]
    }
    myDiagram.model = go.Model.fromJson(full_adder);
    document.getElementById("mySavedModel").value = myDiagram.model.toJson();
}

function load_library() {
    var library = document.getElementById("library_id");
    var selected = library.options[library.selectedIndex].value;
    var upload_library = document.getElementById('upload_library');
    if (selected=='Upload a new library') {
        upload_library.style.display = 'block';
    }
    else {
        upload_library.style.display = 'none';
    }
}

////////////////////////

window.addEventListener('DOMContentLoaded', init);
