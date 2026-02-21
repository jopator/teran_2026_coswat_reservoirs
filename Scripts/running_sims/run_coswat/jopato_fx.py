from datetime import datetime
import shutil
from pathlib import Path
import os
import multiprocessing
import stat

def update_print_prt(file_path, cfg):
    with open(file_path) as f:
        lines = f.readlines()

    out = []
    i = 0
    n = len(lines)

    # 1) General settings
    while i < n:
        line = lines[i]
        stripped = line.strip()

        out.append(line)
        i += 1

        # stop when we hit the object table
        if stripped.startswith("objects"):
            break

        parts = line.split()
        if not parts:
            continue

        # if next line exists, treat this as potential "header row"
        if i < n:
            next_line = lines[i]
            val_parts = next_line.split()

            changed = False
            if val_parts:
                # for each column, if header name is in cfg and is a string, replace its value
                for idx, name in enumerate(parts):
                    if name in cfg and isinstance(cfg[name], str) and idx < len(val_parts):
                        val_parts[idx] = cfg[name]
                        changed = True

            if changed:
                # rebuild the value line simply with spaces
                new_val_line = "\t\t\t".join(val_parts) + "\n"
                out.append(new_val_line)
                i += 1


    # 2) Objects
    col_order = ["daily", "monthly", "yearly", "avann"]

    while i < n:
        line = lines[i]
        parts = line.split()

        if len(parts) < 5:
            out.append(line)
            i += 1
            continue

        name = parts[0]

        if name in cfg and isinstance(cfg[name], dict):
            flags = cfg[name]

            # original flags from file
            orig_flags = parts[1:5]

            new_flags = []
            for col_idx, col in enumerate(col_order):
                if col in flags:
                    new_flags.append(flags[col])
                else:
                    new_flags.append(orig_flags[col_idx])

            # rebuild line with simple spacing; adjust widths if you want nicer alignment
            new_line = f"{name:<28}"
            for v in new_flags:
                new_line += f"{v:<13}"
            out.append(new_line.rstrip() + "\n")
        else:
            out.append(line)

        i += 1

    return out


def copyWeatherFiles(src_dir: Path, dst_dir: Path) -> None:

    cli_files = {
        "wgn": "weather-wgn.cli",
        "wst": "weather-sta.cli",
        "pcp": "pcp.cli",
        "tmp": "tmp.cli",
        "hmd": "hmd.cli",
        "slr": "slr.cli",
        "wnd": "wnd.cli"
    }

    var_files = {
        "pcp": ".pcp",
        "tmp": ".tmp",
        "hmd": ".hmd",
        "slr": ".slr",
        "wnd": ".wnd"
    }

    src = Path(src_dir)
    dst = Path(dst_dir)

    # cli_files
    for fname in cli_files.values():
        fsrc = src / fname
        if fsrc.is_file():
            shutil.copy(fsrc, dst)

    # all files matching the extensions in var_files
    for ext in var_files.values():
        for fsrc in src.glob(f"*{ext}"):
            fdst = dst / fsrc.name
            shutil.copy(fsrc, dst)


def submitCopy(src_dir: str, dst_dir: str):
    copyWeatherFiles(Path(src_dir), Path(dst_dir))
    return src_dir



def sortLexi(region: Path) -> None:
    cli_files = {
        "wgn": "weather-wgn.cli",
        "wst": "weather-sta.cli",
        "pcp": "pcp.cli",
        "tmp": "tmp.cli",
        "hmd": "hmd.cli",
        "slr": "slr.cli",
        "wnd": "wnd.cli"
    }

    for fname in cli_files.values():
        cli_file = region / fname

        if not cli_file.is_file():
            continue

        with open(cli_file) as f:
            lines = f.readlines()

        header = lines[:2]
        data   = lines[2:]

        if fname == "weather-sta.cli":
            data_sorted = sorted(data, key=lambda x: x.split()[0])
            new_header = header  # keep original

        elif fname != "weather-wgn.cli":
            data_sorted = sorted(data)
            new_header = header

        else:
            continue

        with open(cli_file, "w") as f:
            f.writelines(new_header)
            f.writelines(data_sorted)

def submitLexi(dst_dir: str):
    sortLexi(Path(dst_dir))
    return dst_dir


def getBaseName(path_, extension=True):
    if extension:
        fn = os.path.basename(path_)
    else:
        fn = os.path.basename(path_).split(".")[0]
    return(fn)

def runModel(swatDir:str, exePath:str = "/user/brussel/100/vsc10035/VO/data_repo/Executables/swatplus_64", logPath: str | None = None) -> None:
    """
    swatDir: string or pathlike object
    exePath: string or pathlike object
    logPath: string or pathlike object

    Runs the SWAT+ model in parallel
    returns None
    """
    os.chdir(swatDir)

    if logPath is None:
        logPath = f"../../{getBaseName(swatDir, extension=False)}-run.log"

    os.chmod(exePath, os.stat(exePath).st_mode | stat.S_IXUSR)
    os.system(f"{exePath} >> {logPath}")


