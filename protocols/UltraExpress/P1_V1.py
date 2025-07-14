from opentrons import protocol_api
from opentrons import types
from opentrons.protocol_api import SINGLE, ALL
import math
import numpy as np

metadata = {
    "protocolName": "NEBNext UltraExpress® RNA Library Prep Kit_NEB #E3330S/L_Part1_V1",
    "author": "Opentrons",
    "description": "Part 1: mRNA Isolation, Fragmentation and Priming Starting with Total RNA",
}

requirements = {"robotType": "Flex", "apiLevel": "2.20"}


def add_parameters(parameters):
    parameters.add_bool(
        variable_name="DRY_WATER_RUN",
        display_name="Dry or Sample Run",
        description="Do you want to perform a dry run?",
        default=True,
    )
    parameters.add_int(
        variable_name="num_samples",
        display_name="Number of Samples",
        description="The number of samples to process",
        default=94,
        minimum=1,
        maximum=96,
        unit="x Samples",
    )
    parameters.add_str(
        variable_name="p50m_mount",
        display_name="50 µL Pipette Mount",
        description="The mount to attach the 50 µL pipette to",
        choices=[
            {"display_name": "Left", "value": "left"},
            {"display_name": "Right", "value": "right"},
        ],
        default="left",
    )
    parameters.add_str(
        variable_name="p1000m_mount",
        display_name="1000 µL Pipette Mount",
        description="The mount to attach the 1000 µL pipette to",
        choices=[
            {"display_name": "Left", "value": "left"},
            {"display_name": "Right", "value": "right"},
        ],
        default="right",
    )
    parameters.add_float(
        variable_name="flow_rate_aspirate",
        display_name="Default Flow rate aspirate",
        description="Flow rate for the aspirate step",
        default=700,
        minimum=5,
        maximum=716,
        unit="µL/s",
    )
    parameters.add_float(
        variable_name="flow_rate_dispense",
        display_name="Default Flow rate dispense",
        description="Flow rate for the dispense step",
        default=700,
        minimum=5,
        maximum=716,
        unit="µL/s",
    )
    parameters.add_float(
        variable_name="sample_rate",
        display_name="Sample Rate",
        description="Sample mixing rate",
        default=0.6,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="beads_rate",
        display_name="Beads Rate",
        description="Beads mixing rate",
        default=0.2,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="buffer_rate",
        display_name="Buffer Rate",
        description="Buffer mixing rate",
        default=0.8,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="ethanol_rate",
        display_name="Ethanol Rate",
        description="Ethanol mixing rate",
        default=0.8,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="elution_rate",
        display_name="Elution Rate",
        description="Elution mixing rate",
        default=0.4,
        minimum=0,
        maximum=1.5,
    )


def run(protocol: protocol_api.ProtocolContext):
    # ------------------------------------------------------------------ #
    #                             Parameters                             #
    # ------------------------------------------------------------------ #

    # ------------------ Load the run time parameters ------------------ #
    DRY_WATER_RUN = protocol.params.DRY_WATER_RUN
    num_samples = protocol.params.num_samples
    p50m_mount = protocol.params.p50m_mount
    p1000m_mount = protocol.params.p1000m_mount
    flow_rate_aspirate = protocol.params.flow_rate_aspirate
    flow_rate_dispense = protocol.params.flow_rate_dispense
    sample_rate = protocol.params.sample_rate
    beads_rate = protocol.params.beads_rate
    buffer_rate = protocol.params.buffer_rate
    ethanol_rate = protocol.params.ethanol_rate
    elution_rate = protocol.params.elution_rate
    # ------------------------ volume parameters ----------------------- #
    sample_vol = 50
    dT25_beads_vol = 20
    binding_buffer_vol = 50
    tris_and_binding_mm_vol = 100
    wash_buffer_vol = 200
    fragmentation_vol_in = 6.5
    fragmentation_vol_out = 5

    # -------------------- magnetic block parameters ------------------- #
    beads_binding_time = 2

    # ------------------------------------------------------------------ #
    #                           Error Handling                           #
    # ------------------------------------------------------------------ #

    # 1. Check if the number of samples is within the range
    if num_samples > 96 or num_samples < 1:
        raise ValueError("The number of samples should be between 1 and 96")
    # 2. Any of the volumes should be equal to or less than 200 µL
    if (
        dT25_beads_vol > 200
        or binding_buffer_vol > 200
        or tris_and_binding_mm_vol > 200
        or wash_buffer_vol > 200
        or fragmentation_vol_in > 200
        or fragmentation_vol_out > 200
    ):
        raise ValueError("The volume should be equal to or less than 200 µL")

    # ------------------------------------------------------------------ #
    #                             Calculation                            #
    # ------------------------------------------------------------------ #
    # 1. Calculate the column number
    cols = math.ceil(num_samples / 8)
    if num_samples % 8 == 0:
        cols_m = cols
        cols_s = None
    elif num_samples % 8 != 0:
        if cols == 1:
            cols_m = 0
            cols_s = num_samples
        else:
            cols_m = cols - 1
            cols_s = num_samples % 8
    # print(f"cols: {cols}, cols_m: {cols_m}, cols_s: {cols_s}")
    # cols_m and cols_s will be used in the dt25_beads preparation step

    # ------------------------------------------------------------------ #
    #                Tips handling calculation and loading               #
    # ------------------------------------------------------------------ #

    # tip usage calculation
    # Rule 1: Only use 50 µL tips for the fragmentation step.
    # Rule 2: Other steps use 200 µL tips for better beads mixing and removal. It also easier to load/refill tips.

    # --------------------- step by step analysis --------------------- #
    # 1. steps 1A2.1 to 1A2.5: 20 µL beads in, 20 µL supernatant out
    # 1 column for beads addition
    # 1 column reserved for partial column, and 1 column for supernatant removal
    # move plate to a regular slot for the next step
    required_cols_200_1 = 1 + 1 + 1

    # 2. steps 1A2.6 to 1A2.8: 50 µL binding buffer in, 50 µL supernatant out
    # 1 column for binding buffer addition on the top of each well + custom mixing after addition all columns
    # 1 column for supernatant removal
    # move plate to a regular slot for the next step
    required_cols_200_2 = cols + 1

    # 3. steps 1A2.9 to 1A2.14: 50 µL binding buffer in, 100 µL supernatant out after TC
    # cols x columns: Adding buffer + mixing + transferring to tc plate(containing 50 µL sample) in each column
    # TC steps: 2 min at 80°C, 5 min at 25°C, hold at 25°C
    # cols x columns: supernatant removal
    # plate stay on magnetic block for the next step
    required_cols_200_3 = cols + cols

    # 4. steps 1A2.15 to 1A2.17: 200 µL wash buffer in, 200 µL supernatant out
    # cols x columns: Adding wash buffer, no mixing; Or adding on the top of each well, no mixing
    # cols x columns: supernatant removal
    # move plate to a regular slot for the next step
    required_cols_200_4 = 1 + cols  # required_cols_200_4 = cols + cols

    # 5. steps 1A2.18 to 1A2.22: 100 µL Tris buffer+ Binding buffer in, 100 µL supernatant out
    # cols x columns: Adding mixed buffer + mixing
    # TC steps: 2 min at 80°C, 5 min at 25°C, hold at 25°C
    # cols x columns: supernatant removal
    # plate stay on magnetic block for the next step
    required_cols_200_5 = cols + cols

    # 6. steps 1A2.23 to 1A2.25: 200 µL wash buffer in, 200 µL supernatant out
    # cols x columns: Adding wash buffer, no mixing; Or adding on the top of each well, no mixing
    # cols x columns: supernatant removal
    # move plate to a regular slot for the next step
    required_cols_200_6 = 1 + cols  # required_cols_200_6 = cols + cols

    # ------------------- calculate 50 ul tips first ------------------- #

    # 7. steps 1A2.26 to 1A2.30: 6.5 µL fragmentation MM in, 5 µL fragmented RNA out to a new sample_plate
    # cols x columns: Adding fragmentation MM + mixing
    # TC steps: 15 min at 94°C, hold at 4°C
    # cols x columns: fragmented RNA removal (quickly)
    # new sample_plate should be kept on the thermocycler for the next step
    if cols * 2 > 12:
        required_slots_num_50 = 2
        slots_50 = ["A2", "A3"]
        regular_slot_list = ["B2", "B3", "C3"]
    else:
        required_slots_num_50 = 1
        slots_50 = ["A2"]
        regular_slot_list = ["B2", "B3", "C3", "A3"]
    expansion_slots = ["A4", "B4", "C4"]

    # ---------------------- calculate 200 ul tips --------------------- #
    required_cols_200 = (
        required_cols_200_1
        + required_cols_200_2
        + required_cols_200_3
        + required_cols_200_4
        + required_cols_200_5
        + required_cols_200_6
    )
    required_slots_num_200 = math.ceil(required_cols_200 / 12)

    # ------------------- assign slots for 200ul tips ------------------ #
    difference = required_slots_num_200 - len(regular_slot_list)
    if difference <= 0:
        slots_200 = regular_slot_list[:required_slots_num_200]
        slots_200_extra = []
        slots_empty_expansion = []
    else:
        slots_200 = regular_slot_list
        slots_200_extra = expansion_slots[:difference]
        slots_empty_expansion = expansion_slots[difference:]

    print(f"Number of samples: {num_samples} x samples")
    print(f"50ul tips number: {required_slots_num_50} x slots")
    print(f"50ul tips columns: {int(cols*2)}")
    print(f"50ul tips slots: {slots_50}")
    print(f"200ul tips number: {required_slots_num_200} x slots")
    print(f"200ul tips columns: {int(required_cols_200)}")
    print(f"200ul tips slots: {slots_200}")
    print(f"200ul tips extra slots: {slots_200_extra}")
    print(f"Slots not used: {slots_empty_expansion}")

    # ---------------------------- load tips --------------------------- #
    tips50 = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot)
        for slot in slots_50
    ]
    tips200 = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
        for slot in slots_200
    ]
    tips200_extra = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
        for slot in slots_200_extra
    ]

    # ------------------------------------------------------------------ #
    #                               Pipettes                             #
    # ------------------------------------------------------------------ #
    p50m = protocol.load_instrument(
        "flex_8channel_50",
        p50m_mount,
        tips50,
    )
    p1000m = protocol.load_instrument(
        "flex_8channel_1000",
        p1000m_mount,
        tips200,
    )
    p1000m.flow_rate.aspirate = flow_rate_aspirate
    p1000m.flow_rate.dispense = flow_rate_dispense
    # ------------------------------------------------------------------ #
    #                               Modules                              #
    # ------------------------------------------------------------------ #
    thermocycler = protocol.load_module("thermocycler module gen2")
    thermocycler.open_lid()
    thermocycler.set_block_temperature(temperature=4)
    mag = protocol.load_module("magneticBlockV1", "D1")
    temp_mod = protocol.load_module("temperature module gen2", "C1")
    temp_adapter = temp_mod.load_adapter("opentrons_96_deep_well_temp_mod_adapter")
    temp_mod.set_temperature(celsius=4)
    chute = protocol.load_waste_chute()
    # ------------------------------------------------------------------ #
    #                               labware                              #
    # ------------------------------------------------------------------ #
    reagent_plate = temp_adapter.load_labware(
        "nest_96_wellplate_2ml_deep", label="Reagent Plate"
    )
    sample_plate = thermocycler.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "Sample Plate"
    )
    mag_prep_plate = protocol.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "C2", "Mag Prep Plate"
    )
    final_plate_p1 = protocol.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt",
        "D4",
        "Final Plate Part 1",
    )
    waste_res = protocol.load_labware("nest_1_reservoir_195ml", "D2", "Waste Reservoir")
    transition_slot = "C2"
    # ------------------------------------------------------------------ #
    #                           Plate mapping                           #
    # ------------------------------------------------------------------ #

    # _cols represent the 1st well of each column that will be used in the liquid transfer
    # _cols_m represent the 1st well of each column that will be used in the dt25_beads preparation
    # _cols_s represent the leftover wells that will be used in the dt25_beads preparation with partial tip pickup
    # _wells represent all the wells that in mag_prep_plate and sample_plate, used in the liquid labeling
    # _list represent the list of all the wells that will be used in the liquid labelling

    # ------------------------- mag prep plate ------------------------- #
    mag_prep_cols_full = mag_prep_plate.rows()[0][:cols]
    mag_prep_cols_m = mag_prep_plate.rows()[0][:cols_m]
    if cols_s:
        mag_prep_cols_s = mag_prep_plate.columns()[cols_m][:cols_s]
    else:
        mag_prep_cols_s = None
    mag_wells = mag_prep_plate.wells()[:num_samples]

    # print(f"mag_prep_cols_m: {mag_prep_cols_m}")
    # print(f"mag_prep_cols_s: {mag_prep_cols_s}")
    # -------------------------- Sample plate -------------------------- #
    sample_plate_cols = sample_plate.rows()[0][:cols]
    sample_plate_wells = sample_plate.wells()[:num_samples]

    # print(f"sample_plate_cols: {sample_plate_cols}")
    # -------------------------- Reagent Plate ------------------------- #
    dt25_beads = reagent_plate.rows()[0][0]  # A1
    dt25_beads_list = reagent_plate.columns()[0]
    binding_buffer = reagent_plate.rows()[0][1]  # A2, used twice
    binding_buffer_list = reagent_plate.columns()[1]
    washing_buffer_1 = reagent_plate.rows()[0][2:4]  # A3 and A4, reserved two slots
    washing_buffer_1_list = reagent_plate.columns()[2] + reagent_plate.columns()[3]
    tris_and_binding_mm = reagent_plate.rows()[0][4]  # A5
    tris_and_binding_mm_list = reagent_plate.columns()[4]
    washing_buffer_2 = reagent_plate.rows()[0][5:7]  # A6 and A7, reserved two slots
    washing_buffer_2_list = reagent_plate.columns()[5] + reagent_plate.columns()[6]
    fragmentation_mm = reagent_plate.rows()[0][7]  # A8
    fragmentation_mm_list = reagent_plate.columns()[7]
    # ------------------------- waste reservoir ------------------------ #
    waste_well = waste_res.wells()[0]

    # summarize reagent plate map in a list
    reagent_plate_list = [
        "Column 1(A1 of Reagent Plate): dT25 Beads",
        "Column 2(A2 of Reagent Plate): Binding Buffer",
        "Column 3 and 4(A3 and A4 of Reagent Plate): Washing Buffer, 1st round",
        "Column 5(A5 of Reagent Plate): Tris and Binding MM",
        "Column 6 and 7(A6 and A7 of Reagent Plate): Washing Buffer, 2nd round",
        "Column 8(A8 of Reagent Plate): Fragmentation MM",
        "Column 9 to 12: Reserved for Waste",
    ]
    print(f"Reagent Plate List: {reagent_plate_list}")
    # ----------------------- Final Plate Part 1 ----------------------- #
    final_plate_cols = final_plate_p1.rows()[0][:cols]
    final_plate_wells = final_plate_p1.wells()[:num_samples]

    # ------------------------------------------------------------------ #
    #                               liquids                              #
    # ------------------------------------------------------------------ #
    l_locations = [
        dt25_beads_list,
        binding_buffer_list,
        washing_buffer_1_list if cols > 8 else washing_buffer_1_list[:8],
        tris_and_binding_mm_list,
        washing_buffer_2_list if cols > 8 else washing_buffer_2_list[:8],
        fragmentation_mm_list,
    ]
    l_volumes = [
        dT25_beads_vol * cols * 1.2,
        binding_buffer_vol * 2 * cols * 1.2,
        wash_buffer_vol * 8 * 1.2 if cols > 8 else wash_buffer_vol * cols * 1.2,
        tris_and_binding_mm_vol * cols * 1.2,
        wash_buffer_vol * 8 * 1.2 if cols > 8 else wash_buffer_vol * cols * 1.2,
        fragmentation_vol_in * cols * 1.2,
    ]
    liquids = [
        "dT25 Beads",
        "Binding Buffer",
        "Washing Buffer",
        "Tris and Binding MM",
        "Washing Buffer",
        "Fragmentation MM",
    ]
    descriptions = [
        "dT25 Beads",
        "Binding Buffer",
        "Washing Buffer, 1st round",
        "Tris and Binding MM",
        "Washing Buffer, 2nd round",
        "Fragmentation MM",
    ]

    colors_full = [
        "#FF0000",  # Red
        "#0000FF",  # Blue
        "#008000",  # Green
        "#FFFF00",  # Yellow
        "#FFC0CB",  # Pink
        "#800080",  # Purple
        "#FFA500",  # Orange
        "#808080",  # Grey
        "#00FFFF",  # Cyan
        "#FF00FF",  # Magenta
        "#00FF00",  # Lime
        "#000080",  # Navy
        "#800000",  # Maroon
        "#808000",  # Olive
        "#008080",  # Teal
        "#C0C0C0",  # Silver
        "#FF6347",  # Tomato
        "#4682B4",  # SteelBlue
        "#D2691E",  # Chocolate
        "#FF4500",  # OrangeRed
        "#8A2BE2",  # BlueViolet
        "#A52A2A",  # Brown
        "#DEB887",  # BurlyWood
        "#5F9EA0",  # CadetBlue
        "#7FFF00",  # Chartreuse
        "#D2691E",  # Chocolate
        "#FF7F50",  # Coral
        "#6495ED",  # CornflowerBlue
        "#FFF8DC",  # Cornsilk
        "#DC143C",  # Crimson
        "#00FFFF",  # Cyan
        "#00008B",  # DarkBlue
        "#008B8B",  # DarkCyan
        "#B8860B",  # DarkGoldenRod
        "#A9A9A9",  # DarkGray
        "#006400",  # DarkGreen
        "#BDB76B",  # DarkKhaki
        "#8B008B",  # DarkMagenta
        "#556B2F",  # DarkOliveGreen
        "#FF8C00",  # DarkOrange
    ]  # 40 colors

    colors_full = [x.upper() for x in colors_full]

    # make a new color list to match the liquids order.
    # Use the same color for the same liquid.
    # If the liquid is used more than once,
    # use the same color for all instances.
    colors = []
    for i, liquid in enumerate(liquids):
        if liquid not in liquids[:i]:
            colors.append(colors_full[i])
        else:
            colors.append(colors[liquids.index(liquid)])

    for liquid, des, color, v, loc_list in zip(
        liquids, descriptions, colors, l_volumes, l_locations
    ):
        liq = protocol.define_liquid(
            name=str(liquid), description=str(des), display_color=color
        )
        for loc in loc_list:
            loc.load_liquid(liquid=liq, volume=v)

    # label the reserved wells in the mag_prep_plate with 0 volume
    for k, well in enumerate(mag_wells):
        liq = protocol.define_liquid(
            f"dT25 Beads Well x{k+1}",
            "Reserve for dT25 Beads, no liquid when starting",
            display_color=colors_full[len(colors)],
        )
        well.load_liquid(liquid=liq, volume=0)

    # label the wells in the sample_plate with 50 µL volume
    for k, well in enumerate(sample_plate_wells):
        liq = protocol.define_liquid(
            f"RNA sample x{k+1}",
            f"RNA sample x{k+1}, 50 µL",
            display_color=colors_full[len(colors) + 1],
        )
        well.load_liquid(liquid=liq, volume=50)

    # label the reserved wells in the final_plate_p1 with 0 volume
    for k, well in enumerate(final_plate_wells):
        liq = protocol.define_liquid(
            f"Fragmented RNA x{k+1}",
            "Reserve for Fragmented RNA product, Part 1, no liquid when starting",
            display_color=colors_full[len(colors) + 2],
        )
        well.load_liquid(liquid=liq, volume=0)
    # ------------------------------------------------------------------ #
    #                       custom functions, basic                      #
    # ------------------------------------------------------------------ #

    # ------------------------- nonlocal params ------------------------ #
    pause_time = 1

    tip_count_200 = 0
    tiporder_200 = [t1 for i in range(len(tips200)) for t1 in tips200[i].rows()[0]]

    def pause_attention(msg="pause", flash=False):
        """pause the robot, flashes the light, and display a message"""
        nonlocal pause_time
        if flash:
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
        protocol.comment(f"\n\n\nPAUSE X{pause_time}")
        protocol.pause(msg)
        pause_time += 1

    # -------------------------- tip handling -------------------------- #
    def pick_up_200():
        nonlocal tip_count_200, tiporder_200, tips200, tips200_extra
        if tip_count_200 == len(tiporder_200):
            if len(tips200_extra) == 0:
                pause_attention("Replace empty 200ul filter tips")
                tiporder_200.clear()
                for a in slots_200:
                    del protocol.deck[a]
                for c in expansion_slots:
                    del protocol.deck[c]
                tips200 = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
                    for slot in slots_200
                ]
                for a in range(len(tips200)):
                    tiporder_200 += tips200[a].rows()[0]
                tips200_extra = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
                    for slot in slots_200_extra
                ]
            else:  # still have extra 200 ul tips on the expansion slots
                tiporder_200.clear()
                if len(tips200_extra) == 3:
                    protocol.comment(f"\n---move out 1: from {slots_200[0]} to {chute}")
                    protocol.move_labware(tips200[0], chute, use_gripper=True)
                else:
                    protocol.comment(
                        f"\n---move out 1: from {slots_200[0]} to {slots_empty_expansion[0]}"
                    )
                    protocol.move_labware(
                        tips200[0], slots_empty_expansion[0], use_gripper=True
                    )
                for b in range(len(slots_200_extra)):
                    protocol.comment(
                        f"\n---move in {b+1}: from {slots_200_extra[b]} to {slots_200[b]}"
                    )
                    protocol.move_labware(
                        tips200_extra[b], slots_200[b], use_gripper=True
                    )
                    if b + 1 < len(slots_200_extra):
                        protocol.comment(
                            f"---move out {b+2}: from {slots_200[b+1]} to {slots_200_extra[b]}"
                        )
                        protocol.move_labware(
                            tips200[b + 1], slots_200_extra[b], use_gripper=True
                        )
                    tiporder_200 += tips200_extra[b].rows()[0]
                tips200_extra.clear()
            tip_count_200 = 0
        p1000m.pick_up_tip(tiporder_200[tip_count_200])
        tip_count_200 += 1

    def pick_up_50():
        try:
            p50m.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            pause_attention(msg="\n\nReplace empty tipracks before resuming.")
            p50m.reset_tipracks()
            p50m.pick_up_tip()

    def tip_disposal(pip):
        if DRY_WATER_RUN:
            pip.return_tip()
        elif pip.has_tip:
            pip.drop_tip(chute)

    def slow_withdraw(pip, well, z, delay_seconds):
        pip.default_speed /= 40
        if delay_seconds > 0:
            protocol.delay(seconds=delay_seconds)
        pip.move_to(well.top(z))
        pip.default_speed *= 40

    # ------------------------------------------------------------------ #
    #                      custom function, advanced                     #
    # ------------------------------------------------------------------ #

    # -------------------------- custom mixing ------------------------- #
    def custom_mix(pip, mvol, mix_loc, mix_rep, blowout=False, low=False, high=False):
        # Adjustable low and high zoffset for bead mixing
        # The default is false for regular mixing
        if low:
            asp = mix_loc.bottom(low)
        else:
            asp = mix_loc.bottom(1.5)
        if high:
            disp = mix_loc.bottom(high)
        else:
            disp = mix_loc.bottom(2.5)

        # define mixing volume, and use the 2nd smallest value
        if pip == p50m:
            a = 2
        elif pip == p1000m:
            a = 5
        b = 0.8 * pip.tip_racks[0].wells()[0].max_volume
        c = mvol
        numbers = [a, b, c]
        vol = sorted(numbers)[1]

        for _ in range(mix_rep):
            pip.aspirate(vol, asp, rate=sample_rate)
            pip.dispense(vol, disp, push_out=0)
        if blowout:
            pip.flow_rate.blow_out /= 5
            slow_withdraw(pip, mix_loc, z=-3, delay_seconds=1)
            pip.blow_out(mix_loc.top(-3))
            pip.touch_tip(radius=0.75, v_offset=-3, speed=10)
            pip.flow_rate.blow_out *= 5

    # -------------------------- custom delay -------------------------- #
    def custom_delay(name, time):
        if DRY_WATER_RUN:
            for j in np.arange(time, 0, -time):
                msg = (
                    "There are "
                    + str(j)
                    + " seconds left in the "
                    + str(name)
                    + " step"
                )  # noqa
                protocol.delay(seconds=time, msg=msg)
        else:
            for j in np.arange(time, 0, -0.5):
                msg = (
                    "There are "
                    + str(j)
                    + " minutes left in the "
                    + str(name)
                    + " step"
                )  # noqa
                protocol.delay(minutes=0.5, msg=msg)

    # -------------------- well to list distribution ------------------- #
    def well_to_list(
        pip,
        well,
        dest_list,
        transfer_vol,
        liquid_rate,
        z_height,
        pre_mix=False,
        post_mix=False,
        post_mix_rep=0,
        post_mix_vol=0,
        tip_withdrawal=False,
        change_tip=False,
        index=0,
    ):
        for i, d in enumerate(dest_list):
            protocol.comment(
                f"\n---Column x{i+index+1}: {transfer_vol} µL from {well} to {d.parent}"
            )
            if change_tip:
                if pip == p1000m:
                    pick_up_200()
                else:
                    pick_up_50()
            if pre_mix:
                custom_mix(pip, transfer_vol, well, 6)
            if pip == p50m and 1 <= transfer_vol <= 5:
                pip.configure_for_volume(transfer_vol)
            pip.aspirate(transfer_vol, well, rate=liquid_rate)
            pip.dispense(transfer_vol, d.bottom(z_height), rate=liquid_rate)
            if post_mix:
                # set a threshold for the post mix volume in case p50m in a low volume mode
                if pip == p50m and 1 <= transfer_vol <= 5:
                    post_mix_vol = 25 if post_mix_vol > 25 else post_mix_vol
                if post_mix_vol <= 10:
                    custom_mix(
                        pip, post_mix_vol, d, post_mix_rep, blowout=False, low=1, high=1
                    )
                else:
                    custom_mix(pip, post_mix_vol, d, post_mix_rep, blowout=False)
            if tip_withdrawal:
                slow_withdraw(pip, d, z=-3, delay_seconds=0)
            if pip == p50m and 1 <= transfer_vol <= 5:
                pip.configure_for_volume(50)
            if z_height < 0.5 * d.depth:
                pip.blow_out(d.top(-3))
            if change_tip:
                tip_disposal(pip)
            # else:
            #     if i == len(dest_list) - 1:
            #         tip_disposal(pip)

    # --------------------- dt25 beads preparation --------------------- #
    def dt25_beads_prep(change_partial_tip=False):
        nonlocal tip_count_200
        protocol.comment("\n---DT25 Beads Preparation---")
        # 1. Add 20 µL of dT25 beads to each well in the first column
        pick_up_200()
        well_to_list(
            p1000m,
            dt25_beads,
            mag_prep_cols_m,
            dT25_beads_vol,
            beads_rate,
            1.5,
            pre_mix=True,
            post_mix=False,
            post_mix_rep=0,
            post_mix_vol=0,
            tip_withdrawal=True,
            change_tip=False,
        )
        tip_disposal(p1000m)
        if cols_s:
            tips200_available = tips200[0].columns()[1][::-1]
            p1000m.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=tips200)
            for i, (s, d) in enumerate(zip(dt25_beads_list, mag_prep_cols_s)):
                protocol.comment(
                    f"\n---Partial Column x{cols_m+1}, well x{i+1}: {dT25_beads_vol} µL"
                )
                if not p1000m.has_tip:
                    p1000m.pick_up_tip(location=tips200_available[i])
                custom_mix(p1000m, dT25_beads_vol, s, 5)
                p1000m.aspirate(dT25_beads_vol, s, rate=beads_rate)
                p1000m.dispense(dT25_beads_vol, d.bottom(1.5), rate=beads_rate)
                slow_withdraw(p1000m, d, z=-3, delay_seconds=0)
                p1000m.blow_out(d.top(-3))
                if change_partial_tip:
                    if not DRY_WATER_RUN:
                        tip_disposal(p1000m)
                    else:
                        p1000m.drop_tip(chute)
                else:
                    if i == len(mag_prep_cols_s) - 1:
                        if not DRY_WATER_RUN:
                            tip_disposal(p1000m)
                        else:
                            p1000m.drop_tip(chute)
            p1000m.configure_nozzle_layout(style=ALL, tip_racks=tips200)
            tip_count_200 += 1
        protocol.move_labware(mag_prep_plate, mag, use_gripper=True)
        custom_delay("dT25 Beads incubation", time=beads_binding_time)

    # ------------------------- Remove liquids ------------------------- #

    def remove(re_vol, re_s_list, mode="elution", b2slots=False, change_tip=False):

        if mode not in [
            "beads_supernatant",
            "binding_buffer",
            "washing_buffer_1",
            "tris_binding_mm",
            "washing_buffer_2",
            "elution",
        ]:
            raise ValueError("A wrong mode is used")
        if mode == "elution":
            dest_list = final_plate_cols
        else:
            dest_list = [waste_well] * cols
        for i, (s, d) in enumerate(zip(re_s_list, dest_list)):
            protocol.comment(f"\n~~~Column X{i+1}: {mode.title()} Removal")
            protocol.comment(f"---Source: {str(s)}")
            protocol.comment(f"---Destination: {str(d)}")
            if mode == "elution":
                re_rate = elution_rate
                loc = d.bottom(1.5)
            elif mode == "ethanol":
                re_rate = ethanol_rate
                loc = d.top(-3)
            elif mode == "beads_supernatant":
                re_rate = sample_rate
                loc = d.top(-3)
            else:
                re_rate = buffer_rate
                loc = d.bottom(1.5)
            # define volume parameters
            if re_vol > 200:
                raise ValueError("A wrong volume is used")
            elif 150 < re_vol <= 200:
                pip = p1000m
                if not pip.has_tip:
                    pick_up_200()
                pip.aspirate(
                    50,
                    s.bottom().move(types.Point(x=0, y=0, z=0.5 * (s.depth))),
                    rate=0.1,
                )
                slow_withdraw(pip, s, z=3, delay_seconds=1)
                pip.dispense(50, loc, rate=re_rate, push_out=0)
                asp_vol = re_vol - 50
                extra_vol = 20
            elif 20 <= re_vol <= 150:
                pip = p1000m
                asp_vol = re_vol
                extra_vol = 20
            else:
                pip = p50m
                asp_vol = re_vol
                extra_vol = 5
            if asp_vol + extra_vol * 2 > pip.tip_racks[0].wells()[0].max_volume:
                extra_vol = (pip.tip_racks[0].wells()[0].max_volume - asp_vol) / 2
            disp_vol = asp_vol + extra_vol

            if not pip.has_tip:
                if pip == p1000m:
                    pick_up_200()
                else:
                    pick_up_50()
            # air gap for elution
            if mode == "elution":
                pip.move_to(s.top())
                pip.air_gap(extra_vol)
            # regular removal
            pip.aspirate(
                asp_vol,
                s.bottom().move(types.Point(x=0, y=0, z=0.7)),
                rate=0.1,
            )
            # extra removal for buffers and ethanol, air gap built up on the bottom
            if mode != "elution":
                pip.aspirate(
                    extra_vol,
                    s.bottom().move(types.Point(x=0, y=0, z=0.5)),
                    rate=0.05,
                )
            slow_withdraw(pip, s, z=-3, delay_seconds=1)
            pip.dispense(disp_vol, loc, rate=re_rate)  # no push out used
            protocol.delay(seconds=3)
            if mode == "elution":
                slow_withdraw(pip, d, z=-3, delay_seconds=0)
                pip.touch_tip(radius=0.75, v_offset=-3, speed=5)
            else:
                pip.blow_out(loc)
                pip.air_gap(extra_vol)
            if change_tip:
                tip_disposal(pip)
            else:
                if i == len(dest_list) - 1:
                    tip_disposal(pip)
        if b2slots:
            protocol.move_labware(
                re_s_list[0].parent, transition_slot, use_gripper=True
            )
            protocol.comment("\n")

    # ------------------------------ wash ------------------------------ #
    def wash(wash_vol, wash_s_list, addition_change_tip=False, plate_move=False):
        height = sample_plate.wells()[0].depth
        s1 = wash_s_list[0]
        dest_list_1 = sample_plate_cols
        s2 = None
        dest_list_2 = None
        if len(wash_s_list) == 2:
            dest_list_1 = sample_plate_cols[:8]
            s2 = wash_s_list[1]
            dest_list_2 = sample_plate_cols[8:]
        protocol.comment("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        protocol.comment(" Washing Buffer wash starts")
        protocol.comment("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        if not addition_change_tip:
            pick_up_200()
        well_to_list(
            p1000m,
            s1,
            dest_list_1,
            wash_vol,
            ethanol_rate,
            height + 0.5,
            pre_mix=False,
            post_mix=False,
            post_mix_rep=0,
            post_mix_vol=0,
            tip_withdrawal=False,
            change_tip=addition_change_tip,
            index=0,
        )
        if s2 is not None and dest_list_2 is not None:
            well_to_list(
                p1000m,
                s2,
                dest_list_2,
                wash_vol,
                ethanol_rate,
                height + 0.5,
                pre_mix=False,
                post_mix=False,
                post_mix_rep=0,
                post_mix_vol=0,
                tip_withdrawal=False,
                change_tip=addition_change_tip,
                index=8,
            )
        if not addition_change_tip:
            tip_disposal(p1000m)
        # custom_delay("Washing buffer incubation", time=0.5)
        wash_mode = "washing_buffer_2" if plate_move else "washing_buffer_1"
        protocol.comment("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        protocol.comment(f" {wash_mode} Removal starts")
        protocol.comment("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        remove(
            wash_vol,
            sample_plate_cols,
            mode=wash_mode,
            b2slots=plate_move,
            change_tip=True,
        )

    def binding_buffer_addition_to_tc():
        protocol.comment("\n---Binding Buffer Addition---")
        for j, (d1, d2) in enumerate(zip(mag_prep_cols_full, sample_plate_cols)):
            pick_up_200()
            well_to_list(
                p1000m,
                binding_buffer,
                [d1],
                binding_buffer_vol,
                buffer_rate,
                2,
                pre_mix=False,
                post_mix=False,
                post_mix_rep=0,
                tip_withdrawal=False,
                change_tip=False,
                index=j,
            )
            well_to_list(
                p1000m,
                d1,
                [d2],
                binding_buffer_vol,
                buffer_rate,
                2,
                pre_mix=True,
                post_mix=True,
                post_mix_rep=6,
                post_mix_vol=(sample_vol + binding_buffer_vol) * 0.8,
                tip_withdrawal=True,
                change_tip=False,
                index=j,
            )
            tip_disposal(p1000m)
            protocol.comment("\n")
        protocol.comment("\n---Discard the mag_prep_plate to save space---")
        protocol.move_labware(mag_prep_plate, chute, use_gripper=True)

    def tc_steps(
        tc_vol,
        lid_temp,
        t1,
        time1,
        t2,
        time2,
        t3,
        time3,
        h_temp,
        custom_message=False,
    ):
        nonlocal tip_count_200, pause_time
        protocol.comment("\n---Thermocycler incubation---")
        if not DRY_WATER_RUN:
            thermocycler.set_lid_temperature(lid_temp)
        thermocycler.close_lid()
        if not DRY_WATER_RUN:
            thermocycler.set_block_temperature(
                t1, hold_time_minutes=time1, block_max_volume=tc_vol
            )
            if t2 and time2:
                thermocycler.set_block_temperature(
                    t2, hold_time_minutes=time2, block_max_volume=tc_vol
                )
            if t3 and time3:
                thermocycler.set_block_temperature(
                    t3, hold_time_minutes=time3, block_max_volume=tc_vol
                )
            thermocycler.set_block_temperature(h_temp)
            if custom_message and required_slots_num_200 >= 7:
                pause_attention(
                    "Please replace the empty and partial 200ul filter tips before resuming."
                )
                tip_count_200 = 0
        else:
            custom_delay("Dry run delay", time=1)
        thermocycler.open_lid()
        thermocycler.deactivate_lid()
        protocol.move_labware(sample_plate, mag, use_gripper=True)
        custom_delay("after thermocycler incubation", time=beads_binding_time)

    def elution(volume_in, volume_out):
        protocol.comment("\n---Elution Starts---")
        well_to_list(
            p50m,
            fragmentation_mm,
            sample_plate_cols,
            volume_in,
            buffer_rate,
            0.8,
            pre_mix=False,
            post_mix=True,
            post_mix_rep=6,
            post_mix_vol=volume_in * 0.8,
            tip_withdrawal=True,
            change_tip=True,
        )
        protocol.move_labware(sample_plate, thermocycler, use_gripper=True)
        tc_steps(
            tc_vol=volume_in,
            lid_temp=105,
            t1=94,
            time1=15,
            t2=None,
            time2=None,
            t3=None,
            time3=None,
            h_temp=4,
            custom_message=False,
        )
        protocol.move_labware(final_plate_p1, thermocycler, use_gripper=True)
        remove(
            volume_out,
            sample_plate_cols,
            mode="elution",
            b2slots=False,
            change_tip=True,
        )
        thermocycler.close_lid()

    # ------------------------------------------------------------------ #
    #                        protocol starts here                        #
    # ------------------------------------------------------------------ #
    protocol.comment("\n---Protocol starts---\n")
    protocol.comment(
        "\n---1A2.1 to 1A2.5: DT25 Beads Preparation + Remove supernatant---\n"
    )
    dt25_beads_prep(change_partial_tip=True)
    remove(
        dT25_beads_vol,
        mag_prep_cols_full,
        mode="beads_supernatant",
        b2slots=True,
        change_tip=False,
    )

    protocol.comment(
        "\n---1A2.6 to 1A2.8: Binding Buffer Addition(1st round) + Remove supernatant---\n"
    )
    well_to_list(
        p1000m,
        binding_buffer,
        mag_prep_cols_full,
        binding_buffer_vol,
        buffer_rate,
        2,
        pre_mix=False,
        post_mix=True,
        post_mix_rep=6,
        post_mix_vol=binding_buffer_vol * 0.8,
        tip_withdrawal=True,
        change_tip=True,
    )
    protocol.move_labware(mag_prep_plate, mag, use_gripper=True)
    custom_delay("1st Binding Buffer Binding", time=beads_binding_time)
    remove(
        binding_buffer_vol,
        mag_prep_cols_full,
        mode="binding_buffer",
        b2slots=True,
        change_tip=False,
    )

    protocol.comment(
        "\n---1A2.9 to 1A2.14: 50 µL Binding Buffer Addition(2nd) + TC steps+ supernatant removal---\n"
    )
    binding_buffer_addition_to_tc()
    tc_steps(
        tc_vol=binding_buffer_vol + sample_vol,
        lid_temp=90,
        t1=80,
        time1=2,
        t2=25,
        time2=5,
        t3=None,
        time3=None,
        h_temp=25,
        custom_message=True,
    )
    remove(
        binding_buffer_vol + sample_vol,
        sample_plate_cols,
        mode="binding_buffer",
        b2slots=False,
        change_tip=True,
    )  # binding buffer + sample supernatant removal
    protocol.comment("\n---1A2.15 to 1A2.17: 1st 200 µL Wash Buffer---\n")
    wash(wash_buffer_vol, washing_buffer_1, addition_change_tip=False)
    protocol.move_labware(sample_plate, transition_slot, use_gripper=True)

    protocol.comment(
        "\n---1A2.18 to 1A2.22: tris and binding MM addition + TC steps+ supernatant removal---\n"
    )
    protocol.comment("\n---Tris and Binding MM Addition---\n")
    well_to_list(
        p1000m,
        tris_and_binding_mm,
        sample_plate_cols,
        tris_and_binding_mm_vol,
        buffer_rate,
        2,
        pre_mix=False,
        post_mix=True,
        post_mix_rep=6,
        post_mix_vol=tris_and_binding_mm_vol * 0.8,
        tip_withdrawal=False,
        change_tip=True,
    )
    protocol.move_labware(sample_plate, thermocycler, use_gripper=True)
    tc_steps(
        tc_vol=tris_and_binding_mm_vol,
        lid_temp=90,
        t1=80,
        time1=2,
        t2=25,
        time2=5,
        t3=None,
        time3=None,
        h_temp=25,
        custom_message=False,
    )
    thermocycler.set_block_temperature(4)  # prepare for final elution
    remove(
        tris_and_binding_mm_vol,
        sample_plate_cols,
        mode="tris_binding_mm",
        b2slots=False,
        change_tip=True,
    )
    protocol.comment("\n---1A2.23 to 1A2.25: 2nd 200 µL Wash Buffer---\n")
    wash(
        wash_buffer_vol,
        washing_buffer_2,
        addition_change_tip=False,
        plate_move=True,
    )

    protocol.comment(
        "\n---1A2.26 to 1A2.30: 6.5 µL fragmentation MM in, 5 µL fragmented RNA out to a new sample_plate---\n"
    )
    elution(fragmentation_vol_in, fragmentation_vol_out)
    protocol.move_labware(sample_plate, "D4", use_gripper=True)
    pause_attention(
        "\n---The Part 1 protocol has been completed. Press 'Resume' button to finish the protocol---\n"
        "Prepare the reagents for Part 2 immediately after the Part 1 protocol is completed.\n"
        "TC and temp_mode will stay on 4 degrees"
    )

    # thermocycler.open_lid()
    # protocol.move_labware(final_plate_p1, transition_slot, use_gripper=True)
    protocol.comment("\n---Protocol ends---\n")

    # protocol.comment(f"final_tip_count_200: {tip_count_200}")
    # protocol.comment(f"pause_time: {pause_time-1}")
