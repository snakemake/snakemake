'use strict';

class MetaData extends React.Component {
    constructor(props) {
        super(props);
    }

    //TODO: parse recursively here (with inner render function)

    render() {
        let metadatalist = [];
        for (const [key, value] of Object.entries(metadata)) {

            metadatalist.push(...[
                e(
                    ListHeading,
                    { key: `${key}-heading`, text: key }
                ),
                e(
                    ListItem,
                    {
                        key: `${key}`,
                        className: "p-1",
                        dangerouslySetInnerHTML: { __html: value }
                    }
                )
            ]);
        }
        console.log(metadatalist);
        return e(
            "ul",
            { },
            metadatalist
        )
    }
}
