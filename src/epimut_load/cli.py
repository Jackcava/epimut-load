import argparse

# from epimut_load.core import run_pipeline
from epimut_load.core import run_pipeline


def main():

    # parser = argparse.ArgumentParser()

    # parser.add_argument(
    #     "--config",
    #     required=True
    # )

    # args = parser.parse_args()

    # run_pipeline(args.config)


    parser = argparse.ArgumentParser(
        prog="epimut-load"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    # run
    run_parser = subparsers.add_parser(
        "run",
        help="Run EML analysis"
    )

    run_parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file"
    )

    args = parser.parse_args()

    if args.command == "run":
        run_pipeline(args.config)



if __name__ == "__main__":
	main()