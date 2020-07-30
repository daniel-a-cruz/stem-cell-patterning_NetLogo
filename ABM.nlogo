extensions [vid csv]
breed [ cells cell ]

; Global counters for gains and loss of cells at specific states over the course of steps (aka ticks)
globals [ altcount bhGain bhLoss blGain blLoss nhGain nhLoss ghGain ghLoss ]

turtles-own [ motion pluri divTrack diffTrack update FGFR ERK GATA6 NANOG ]
patches-own [ FGF4 ]

; Labeled "Initialize" on interface
to setup

  ; Resets steps and clears all turtles, patches, plots, etc.
  clear-all

  ; Resets recorder for visualization output
  vid:reset-recorder

  ; Counter which tracks how many cells have equal levels of NANOG and GATA6 (i.e. both high or both low)
  set altcount 0

  ; Counters which track the gains/losses of cells that have both NANOG and GATA6 set to high (1)
  set bhGain 0
  set bhLoss 0

  ; Counters which track the gains/losses of cells that have both NANOG and GATA6 set to low (0)
  set blGain 0
  set blLoss 0

  ; Counters which track the gains/losses of cells that have just NANOG set to high (1)
  set nhGain 0
  set nhLoss 0

  ; Counters which track the gains/losses of cells that have just GATA6 set to high (1)
  set ghGain 0
  set ghLoss 0

  ; The shape of all cells is a circle
  set-default-shape turtles "circle"

  ; Patches have their initial value of FGF4 set to a random value up to the user-defined "maxFGF4"
  ask patches [ set FGF4 random maxFGF4 ]

  ; Creates desired amount of pluripotent cells
  ; This is the sum of NANOG High and GATA6 cells, set from the interface
  create-cells (numNanogHigh + numGataHigh) [

    ;
    set update random pathway-update

    ; Each cell is placed at a different spot in the environment
    setxy random-xcor random-ycor

    ; motion is used for conditionals to tell if a cell is able to move
    set motion true

    ; pluri is used for conditionals to tell if a cell is pluripotent
    set pluri true

    ; divTrack is used to determine when a cell divides
    ; A random value is assigned between 0 and interface threshold
    set divTrack random pluriMitosisThreshold

    ; diffTrack is used to determine when a cell differentiates
    ; A random value is assigned between 0 and interface threshold
    set diffTrack random diffThresh

    ; Cells have their initial values of FGFR and ERK set to a random value
    set FGFR random 2
    set ERK random 2

    ; The following determines that a cell has a low (0) value of GATA6 and a high (1) value for NANOG
    set GATA6 0
    set NANOG 1

    ; Lime (~green) is used to mark NANOG High pluripotent cells
    set color lime
  ]

  ; From the cells created above, a random "numGataHigh" number of them are converted to GATA6 High
  ask n-of numGataHigh cells [

    ; The following determines that a cell has a high (1) value of GATA6 and a low (0) value for NANOG
    set GATA6 1
    set NANOG 0

    ; White is used to mark GATA6 High pluripotent cells
    set color white
  ]

  ; Resets the number of steps in preparation for the simulation to begin
  reset-ticks
end

; Labeled "Step" on interface
to go

  ; Resets all counters below to zero in preparation for recording new changes
  set bhGain 0
  set bhLoss 0
  set blGain 0
  set blLoss 0
  set nhGain 0
  set nhLoss 0
  set ghGain 0
  set ghLoss 0

  ask cells [

    ;
    set update (update + 1) mod (pathway-update + 1)

    ; Increases divTrack value by 5 to represent that a cell is closer to division
    set divTrack divTrack + 5

    ; Calls the single-cell-move function, which controls how cells move around the environment
    single-cell-move

    ; Calls the function FE-pathway, which simulates the FGF/ERK pathway for a cell
    if update = pathway-update [ FE-pathway ]

    ; Calls the diff-low-surrounded function, which simulates differentiation caused by being surrounded by differentiated cells
    diff-low-surrounded

    ; Calls the cell-division function, which simulates mitosis for pluripotent and differentiated cells
    cell-division

    ; Calls the color-update function, which updates the color of a cell based on its state
    color-update

  ]

  ; Diffuses 50% of the amount of FGF4 on each patch to its eight neighboring patches
  diffuse FGF4 0.5

  tick

  ; Yellow is used to mark pluripotent cells with both GATA6 and NANOG high
  let ylw count cells with [ color = yellow ]

  ; Sky (~blue) is used to mark pluripotent cells with both GATA6 and NANOG low
  let blu count cells with [ color = sky ]

  set altcount (ylw + blu)
end

; Cone movement, currently not in use
to cone-move [ dist ]
  if dist >= .25 [
    let obs count other cells in-cone 1 60
    ifelse obs < 1 [
      fd dist
    ]
    [
      let angle random-float 30
      let temp random 2
      ifelse temp = 1 [
        rt 30 + angle
      ]
      [
        lt 30 + angle
      ]
      move dist / 2
    ]
  ]
end

;
to move [ dist ]
  let iter 10

  while [iter > 0] [
    fd dist / 10
    let obs count other cells in-cone 1 60
    ifelse obs >= 1 [
      fd -1 * (dist / 10)
      set iter 0
    ]
    [
    set iter (iter - 1)
    ]
  ]
end

; This function controls the movement of an individual cell over the course of one step
to single-cell-move
  if motion [

    ; Differentiated cells move randomly
    ifelse not pluri [
      rt random-float 360
      move cellSpeed

      ; If a cell is differentiated and is within a radius of 1 from another differentiated cell, it will stop moving
      if any? other cells with [ pluri = false ] in-radius 1 [
        set motion false
      ]
    ]
    ; Below, we handle the cases when a cell is pluripotent
    [
      ; If the cell is GATA6 High (and NANOG Low)
      ifelse (GATA6 = 1) and (NANOG = 0) [

        ; Guye Model Movement
        ifelse guyeMove [

          ; Pluirpotent GATA6 High cells localize to the nearest differentiated cell if any exist
          let nearestDiff min-one-of cells in-radius inter-distance with [ pluri = false ] [ distance myself ]

          ifelse nearestDiff != nobody [
            face nearestDiff
            move cellSpeed
          ]
          ; Otherwise, they move at random
          [
            rt random-float 360
            move cellSpeed
          ]
        ]
        ; If Guye model movement is not active, pluripotent GATA6 High cells will just move at random
        [
          rt random-float 360
          move cellSpeed
        ]
      ]
      [
        ; If the cell is NANOG High (and GATA6 Low) and NANOG High cells cluster
        ifelse (GATA6 = 0) and (NANOG = 1) [

          ; NANOG High cells cluster
          ifelse nClustering [

            ; Pluripotent NANOG High cells localize to each other based on nearest neighbors (if any exist)
            let nearestLow min-one-of cells in-radius inter-distance with [ NANOG = 0 ] [ distance myself ]
            ifelse nearestLow != nobody [
              face nearestLow
              move cellSpeed
            ]
            ; Otherwise, they move at random
            [
              rt random-float 360
              move cellSpeed
            ]
          ]
          ; If NANOG High clustering is not active, pluripotent NANOG High cells will just move at random
          ; This case also causes all other cells to move at random
          [
            rt random-float 360
            move cellSpeed
          ]

          ; pluripotent NANOG High cells stop moving if they are near other pluripotent NANOG High cells
          if any? other cells with [ (NANOG = 1) and (GATA6 = 0) and (pluri = true) ] in-radius 1 [
              set motion false
          ]
        ]
        ; All other cells to move at random
        [
          rt random-float 360
          move cellSpeed
        ]
      ]
    ]
  ]
end

; This function changes a cell's variables to represent that it has differentiated
to differentiate
  ; Cells may have stopped moving (i.e clustered) before differentiating, so their motion is reset to true
  set motion true

  set pluri false
  set color red
  set GATA6 1
  set NANOG 0
end


; This function simulations the FGF/ERK pathway and its effects on the cells within the model.
; It also updates the global counters (per time step) based on the state changes that each cell undergoes
to FE-pathway

  ; First, temporary variables are used to store the current values of a cell and the FGF4 value from the patch on which it is located
  let boolFGF4 0
  let tempFGFR FGFR
  let tempERK ERK
  let tempNANOG NANOG
  let tempGATA6 GATA6

  ; This sub-function determines whether the patch on which a cell is located has at least 1 "unit" of FGF4
  ask patch-here [
    if FGF4 >= 1 [
      ; If the patch does, the value of FGF4 will count as "High" or 1
      ; Otherwise, it is set to "Low" or 0 by default
      set boolFGF4 1
    ]
    ; Boolean function for FGF4: If a cell is NANOG High (regardless of other values), then it produces FGF4 onto the patch on which it is located
    if tempNANOG = 1 [
      set FGF4 (FGF4 + 1)
    ]
  ]

  ; Next, the remaining Boolean functions are updated based on their current values
  set FGFR (boolFGF4 * tempGATA6)
  set ERK tempFGFR
  set GATA6 (1 + tempNANOG + tempNANOG * tempGATA6) mod 2
  set NANOG (1 + tempERK + tempGATA6 + tempERK * tempGATA6) mod 2

  ; This conditional represents the "expenditure" of FGF4 to when it is received by FGFR
  if (tempFGFR = 0) and (FGFR = 1) [
    ask patch-here [ set FGF4 FGF4 - 1 ]
  ]

  ifelse (GATA6 = 1) and pluri [

    ifelse tempGATA6 = 1 [

      set diffTrack diffTrack + 1
      if diffTrack >= diffThresh [ differentiate ]
    ]
    [
      set diffTrack random diffThresh
    ]
  ]
  [
    set diffTrack 0
  ]

  ; This conditional was added to prevent cells that are not NANOG High and GATA6 Low from clustering if there was a change in these values
  if tempNANOG != NANOG and pluri [
    set motion true
  ]

  if stocNANOG [

    let stoch random 2

    if (GATA6 = 0) and pluri and (stoch = 1) [
      set NANOG (NANOG + 1) mod 2
    ]
  ]

  ; The conditionals after this point are used to update the global counters which track each kind of cell. These counters were introduced above

  ; NANOG High tracking

  if (NANOG = 1) and (GATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set nhGain nhGain + 1
  ]

  if (tempNANOG = 1) and (tempGATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set nhLoss nhLoss + 1
  ]

  ; GATA6 High tracking

  if (NANOG = 0) and (GATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set ghGain ghGain + 1
  ]

  if (tempNANOG = 0) and (tempGATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6)) [
    set ghLoss ghLoss + 1
  ]

  ; Both High tracking

  if (NANOG = 1) and (GATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set bhGain bhGain + 1
  ]

  if (tempNANOG = 1) and (tempGATA6 = 1) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set bhLoss bhLoss + 1
  ]

  ; Both Low tracking

  if (NANOG = 0) and (GATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set blGain blGain + 1
  ]

  if (tempNANOG = 0) and (tempGATA6 = 0) and ((tempNANOG != NANOG) or (tempGATA6 != GATA6))
  [
    set blLoss blLoss + 1
  ]

end

; Function which allows pluripotent GATA6 low cells to differentiate if they are "surrounded" by differentiated cells (and if the diffLowSurr is enabled)
; Differentiation will be implemented by converting a cell into being GATA6 High
to diff-low-surrounded
  if diffLowSurr [

    if (GATA6 = 0) and pluri [

      ; The function checks the neighbors within a radius of 1 to see if there are enough differentiated neighbors near this cell
      let crowd count other cells with [ pluri = false ] in-radius 1

      ; Maximum number of cells at a distance of 1 "unit" is 6 (without overlap)
      if crowd >= 4 [

        ; If so, then the cell becomes GATA6 High
        set NANOG 0
        set GATA6 1
        set diffTrack random diffThresh
        set motion true

        ifelse NANOG = 0 [ set blLoss blLoss + 1 ] [ set nhLoss nhLoss + 1 ]

        set ghGain ghGain + 1
      ]
    ]
  ]
end

; Function for cell division for pluripotent NANOG High and differentiated cells
to cell-division

  let crowd count other cells in-radius 1

  ; Pluripotent NANOG High cells and differentiated cells can divide if they are clustered with similar cells (and thus have stopped moving).
  ; However, if a cell is completely surrounded, then it cannot divide. The maximum number of cells directly next to a cell is 6 (without overlap).
  if (crowd < 6) and (not motion) [

    ; If a cell is differentiated, the differentiated mitosis threshold is checked against to see if the cell should divide
    ; Otherwise, the pluripotent mitosis threshold is used instead.
    if (not pluri and (divTrack > diffMitosisThreshold)) or (pluri and (divTrack > pluriMitosisThreshold)) [
      set divTrack (divTrack / 2)

      ; The new cell appears on top of the original cell
      hatch 1 [

        set diffTrack random diffThresh
        let divLoop 12
        let preXcor xcor
        let preYcor ycor
        rt random-float 360

       while [ divLoop > 0 ] [

          let obs count other cells in-cone 1 60

          ifelse (obs <= 1) [
            move cellSpeed
            set divLoop 0
          ]
          [
            rt 30
            set divLoop (divLoop - 1)
          ]
        ]

        if (preXcor = xcor) and (preYcor = ycor) [ die ]
      ]
    ]
  ]
end

; Function to update a cell's color based on its GATA6 and NANOG levels, and whether or not it's differentiated
to color-update
  ifelse not pluri [
    ; Differentiated cells are red
    set color red
  ]
  ; The following are all colors for pluripotent cells
  [
    ; NANOG High is lime (~green)
    if (NANOG = 1) and (GATA6 = 0) [
      set color lime
    ]
    ; GATA6 High is white
    if (NANOG = 0) and (GATA6 = 1) [
      set color white
    ]
    ; Both High is yellow
    if (NANOG = 1) and (GATA6 = 1) [
      set color yellow
    ]
    ; Both low is sky (~blue)
    if (NANOG = 0) and (GATA6 = 0) [
      set color sky
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
275
10
889
625
-1
-1
6.0
1
15
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
steps
2.0

SLIDER
5
130
265
163
numGataHigh
numGataHigh
0
1000
0.0
10
1
NIL
HORIZONTAL

SLIDER
5
200
265
233
numNanogHigh
numNanogHigh
0
1000
1000.0
10
1
NIL
HORIZONTAL

BUTTON
900
430
973
463
Initialize
setup\n\nif record-type = \"View\" [\nvid:start-recorder\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:start-recorder\nvid:record-interface\n]    
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1070
430
1155
463
Run
go\n\nif record-type = \"View\" [\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:record-interface\n] 
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

MONITOR
700
635
780
680
Differentiated
count cells with [ color = red ]
17
1
11

SLIDER
5
270
265
303
pluriMitosisThreshold
pluriMitosisThreshold
0
500
75.0
5
1
NIL
HORIZONTAL

SLIDER
5
340
265
373
diffMitosisThreshold
diffMitosisThreshold
0
500
100.0
5
1
NIL
HORIZONTAL

MONITOR
280
635
360
680
NANOG High
count cells with [ color = lime ]
17
1
11

BUTTON
982
430
1062
463
Step
go\n\nif record-type = \"View\" [\nvid:record-view\n]\n\nif record-type = \"Interface\" [\nvid:record-interface\n] 
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1045
475
1155
515
Save Video
if record-type = \"View\" or record-type = \"Interface\" [\nvid:save-recording user-new-file\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
900
475
1035
520
record-type
record-type
"None" "View" "Interface"
2

MONITOR
385
635
465
680
GATA6 High
count cells with [ color = white ]
17
1
11

TEXTBOX
75
10
180
40
Parameters
20
0.0
1

SLIDER
5
60
265
93
cellSpeed
cellSpeed
0
5
1.0
.5
1
NIL
HORIZONTAL

TEXTBOX
965
10
1080
40
Assumptions
20
0.0
1

SWITCH
900
70
1155
103
guyeMove
guyeMove
0
1
-1000

TEXTBOX
905
45
1150
63
GATA6 Movement: Random or Guye?
13
0.0
1

TEXTBOX
905
200
1150
231
Pluripotent GATA6 Low cells differentiate within differentiated subpopulation?
13
0.0
1

SWITCH
900
240
1155
273
diffLowSurr
diffLowSurr
0
1
-1000

TEXTBOX
985
390
1060
410
Controls
20
0.0
1

TEXTBOX
10
40
260
58
Cell speed
13
0.0
1

TEXTBOX
10
110
260
128
Number of initial GATA6 High Cells
13
0.0
1

TEXTBOX
10
180
260
200
Number of initial NANOG High Cells
13
0.0
1

TEXTBOX
10
320
260
338
Mitosis threshold for pluripotent cells
13
0.0
1

TEXTBOX
10
250
260
268
Mitosis threshold for differentiated cells\t
13
0.0
1

TEXTBOX
10
390
250
421
Threshold (in steps) for differentiating cells
13
0.0
1

SLIDER
5
410
265
443
diffThresh
diffThresh
0
45
36.0
3
1
NIL
HORIZONTAL

MONITOR
490
635
570
680
Both High
count cells with [ color = yellow ]
17
1
11

TEXTBOX
1265
10
1345
31
Monitors
20
0.0
1

MONITOR
595
635
677
680
Both Low
count cells with [ color = sky ]
17
1
11

TEXTBOX
285
685
355
705
Green
15
0.0
1

TEXTBOX
390
685
460
705
White
15
0.0
1

TEXTBOX
495
685
565
705
Yellow
15
0.0
1

TEXTBOX
600
685
670
703
Blue
15
0.0
1

TEXTBOX
705
685
775
703
Red
15
0.0
1

MONITOR
1315
115
1430
160
GATA6 High Loss
ghLoss
17
1
11

MONITOR
1315
60
1430
105
NANOG High Loss
nhLoss
17
1
11

MONITOR
1175
60
1295
105
NANOG High Gain
nhGain
17
1
11

MONITOR
1175
115
1295
160
GATA6 High Gain
ghGain
17
1
11

MONITOR
805
635
885
680
Alt. Cells
altcount
17
1
11

TEXTBOX
810
685
880
721
Blue + Yellow
15
0.0
1

MONITOR
1175
170
1295
215
Both High Gain
bhGain
17
1
11

MONITOR
1315
170
1430
215
Both High Loss
bhLoss
17
1
11

MONITOR
1175
225
1295
270
Both Low Gain
blGain
17
1
11

MONITOR
1315
225
1430
270
Both Low Loss
blLoss
17
1
11

TEXTBOX
905
120
1140
138
NANOG High clustering?
13
0.0
1

SWITCH
900
145
1155
178
nClustering
nClustering
0
1
-1000

SLIDER
5
485
265
518
maxFGF4
maxFGF4
1
5
5.0
1
1
NIL
HORIZONTAL

TEXTBOX
10
465
260
483
Maximum initial FGF4 amount on a patch
13
0.0
1

SLIDER
5
555
265
588
pathway-update
pathway-update
0
5
0.0
1
1
NIL
HORIZONTAL

TEXTBOX
10
535
255
553
Number of steps for pathway update
13
0.0
1

SLIDER
5
630
265
663
inter-distance
inter-distance
0
10
5.0
1
1
NIL
HORIZONTAL

TEXTBOX
15
605
260
623
Interaction distance for movement
13
0.0
1

SWITCH
900
325
1155
358
stocNANOG
stocNANOG
1
1
-1000

TEXTBOX
905
300
1150
318
Is NANOG stochastic?
13
0.0
1

CHOOSER
900
535
1035
580
export-cells
export-cells
"GATA6 High" "NANOG High" "Other" "All"
3

BUTTON
1045
535
1155
575
Save CSV
if export-cells = \"GATA6 High\" [\n  csv:to-file user-new-file [ (list xcor ycor) ] of cells with [ color = white or color = red ]\n]\n\nif export-cells = \"NANOG High\" [\n  csv:to-file user-new-file [ (list xcor ycor) ] of cells with [ color = lime ]\n]\n\nif export-cells = \"Other\" [\n  csv:to-file user-new-file [ (list xcor ycor) ] of cells with [ color = sky or color = yellow ]\n]\n\nif export-cells = \"All\" [\n  csv:to-file user-new-file [ (list xcor ycor) ] of cells\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
